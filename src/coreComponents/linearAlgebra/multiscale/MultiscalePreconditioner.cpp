/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiscalePreconditioner.cpp
 */

#include "MultiscalePreconditioner.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/LevelBuilderBase.hpp"

namespace geosx
{

template< typename LAI >
MultiscalePreconditioner< LAI >::MultiscalePreconditioner( LinearSolverParameters params,
                                                           DomainPartition & domain )
  : Base(),
  m_params( std::move( params ) ),
  m_domain( domain )
{}

template< typename LAI >
MultiscalePreconditioner< LAI >::~MultiscalePreconditioner() = default;

template< typename LAI >
void MultiscalePreconditioner< LAI >::printLevelInfo() const
{
  constexpr char const lineFormat[] = "{:>2}  {:>10}  {:>12}  {:>7.4f}  {:>7.2f}  {:>7.2f}\n";
  constexpr char const headFormat[] = "{:>2}  {:>10}  {:>12}  {:>7}  {:>7}  {:>7}\n";

  std::ostringstream os;
  string const header = GEOSX_FMT( headFormat, "L", "rows", "entries", "sparse", "nnz/row", "ratio" );
  os << "\nOperators:\n" << header;
  os << string( header.length() - 1, '=' ) << "\n";

  globalIndex totalNumRows = 0;
  globalIndex totalNumNonzeros = 0;
  globalIndex totalMemory = 0;

  for( size_t levelIndex = 0; levelIndex < m_levels.size(); ++levelIndex )
  {
    Level const & level = m_levels[levelIndex];
    globalIndex const nrow = level.matrix->numGlobalRows();
    globalIndex const nnz = level.matrix->numGlobalNonzeros();
    globalIndex const prevNrow = levelIndex > 0 ? m_levels[levelIndex-1].matrix->numGlobalRows() : nrow;
    os << GEOSX_FMT( lineFormat, levelIndex, nrow, nnz, real64( nnz ) / ( nrow * nrow ), real64( nnz ) / nrow, real64( prevNrow ) / nrow );

    totalNumRows += nrow;
    totalNumNonzeros += nnz;
    totalMemory += nnz;
    if( levelIndex > 0 )
    {
      totalMemory += level.builder->prolongation().numGlobalNonzeros() + level.builder->restriction().numGlobalNonzeros();
    }
  }

  Matrix const & fineMat = *m_levels[0].matrix;
  constexpr char const compFormat[] = "  {:>8} = {:>6.4f}\n";
  os << "\nComplexities:\n";
  os << GEOSX_FMT( compFormat, "grid", real64( totalNumRows ) / fineMat.numGlobalRows() );
  os << GEOSX_FMT( compFormat, "operator", real64( totalNumNonzeros ) / fineMat.numGlobalNonzeros() );
  os << GEOSX_FMT( compFormat, "memory", real64( totalMemory ) / fineMat.numGlobalNonzeros() );

  GEOSX_LOG_RANK_0( os.str() );
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::createLevels( Matrix const & mat,
                                                    DofManager const & dofManager )
{
  GEOSX_MARK_FUNCTION;
  m_levels.clear();

  // create fine level
  {
    m_levels.emplace_back();
    Level & fine = m_levels[0];
    fine.builder = multiscale::LevelBuilderBase< LAI >::createInstance( m_params.multiscale.fieldName + "_L0", m_params.multiscale );
    fine.builder->initializeFineLevel( m_domain, dofManager, m_params.multiscale.fieldName, mat.comm() );
    fine.matrix = &mat;
  }

  // create coarse levels
  LinearSolverParameters::Multiscale params = m_params.multiscale;
  for( integer levelIndex = 1; levelIndex < m_params.multiscale.maxLevels; ++levelIndex )
  {
    m_levels.emplace_back();
    Level & coarse = m_levels[levelIndex];
    string const name = GEOSX_FMT( "{}_L{}", m_params.multiscale.fieldName, levelIndex );
    coarse.builder = multiscale::LevelBuilderBase< LAI >::createInstance( name, m_params.multiscale );
    coarse.builder->initializeCoarseLevel( *m_levels[levelIndex - 1].builder );
    coarse.matrix = &coarse.builder->matrix();

    if( coarse.matrix->numGlobalRows() <= m_params.multiscale.coarsening.minGlobalDof )
    {
      // Prevent further coarsening globally by truncating level hierarchy
      break;
    }
    if( coarse.matrix->numLocalRows() <= m_params.multiscale.coarsening.minLocalDof )
    {
      // Prevent further coarsening locally by setting ratio to 1
      params.coarsening.ratio.setValues< serialPolicy >( 1.0 );
    }
  }

  // create smoothers
  for( size_t levelIndex = 0; levelIndex < m_levels.size() - 1; ++levelIndex )
  {
    Level & level = m_levels[levelIndex];
    // TODO: smoother options from input
    LinearSolverParameters smoother_params;
    smoother_params.preconditionerType = m_params.multiscale.smootherType;
    if( m_params.multiscale.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::pre
        || m_params.multiscale.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::both )
    {
      level.presmoother = LAI::createPreconditioner( smoother_params );
    }
    if( m_params.multiscale.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::post
        || m_params.multiscale.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::both )
    {
      level.postsmoother = LAI::createPreconditioner( smoother_params );
    }
    // TODO: pre/post smoother could be the same object
  }

  // create vectors
  for( Level & level : m_levels )
  {
    level.rhs.create( level.matrix->numLocalRows(), level.matrix->comm() );
    level.sol.create( level.matrix->numLocalRows(), level.matrix->comm() );
    level.tmp.create( level.matrix->numLocalRows(), level.matrix->comm() );
  }

  // create coarse solver
  {
    LinearSolverParameters coarseParams;
    coarseParams.preconditionerType = m_params.multiscale.coarseType;
    m_coarse_solver = LAI::createPreconditioner( coarseParams );
  }
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::setup( Matrix const & mat )
{
  GEOSX_MARK_FUNCTION;

  Base::setup( mat );

  auto const logMessage = [&]( integer const minLevel, string const & msg )
  {
    GEOSX_LOG_RANK_0_IF( m_params.logLevel >= minLevel,
                         "[Multiscale] " << m_params.multiscale.fieldName << ": " << msg );
  };

  if( !m_initialized )
  {
    GEOSX_ERROR_IF( mat.dofManager() == nullptr,
                    GEOSX_FMT( "[Multiscale] {}: DofManager is required", m_params.multiscale.fieldName ) );
    logMessage( 3, "creating level hierarchy" );
    createLevels( mat, *mat.dofManager() );
    m_initialized = true;
  }

  // compute level operators
  for( size_t levelIndex = 1; levelIndex < m_levels.size(); ++levelIndex )
  {
    logMessage( 3, GEOSX_FMT( "computing operator for coarse level {}", levelIndex ) );
    m_levels[levelIndex].builder->compute( *m_levels[levelIndex - 1].matrix );
  }

  // compute level smoothers
  logMessage( 3, "setting up smoothers" );
  for( size_t levelIndex = 0; levelIndex < m_levels.size() - 1; ++levelIndex )
  {
    Level & level = m_levels[levelIndex];
    if( level.presmoother )
    {
      level.presmoother->setup( *level.matrix );
    }
    if( level.postsmoother )
    {
      level.postsmoother->setup( *level.matrix );
    }
  }

  // setup coarse solver
  logMessage( 3, "setting up coarse solver" );
  m_coarse_solver->setup( *m_levels.back().matrix );
  logMessage( 3, "setup completed" );

  logMessage( 1, "statistics:" );
  if( m_params.logLevel >= 1 )
  {
    printLevelInfo();
  }
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::apply( Vector const & src,
                                             Vector & dst ) const
{
  // TODO: remove hardcoded V-cycle, abstract into a separate component

  m_levels[0].rhs.copy( src );
  int const numLevels = LvArray::integerConversion< int >( m_levels.size() );

  // down phase
  for( int levelIndex = 0; levelIndex < numLevels - 1; ++levelIndex )
  {
    Level const & fine = m_levels[levelIndex];
    Level const & coarse = m_levels[levelIndex + 1];
    fine.sol.zero();
    if( fine.presmoother )
    {
      for( integer s = 0; s < m_params.multiscale.numSmootherSweeps; ++s )
      {
        fine.presmoother->apply( fine.rhs, fine.tmp );
        fine.sol.axpy( 1.0, fine.tmp );
        fine.matrix->residual( fine.tmp, fine.rhs, fine.rhs );
      }
    }
    coarse.builder->restriction().apply( fine.rhs, coarse.rhs );
  }

  // coarse level solve
  m_coarse_solver->apply( m_levels.back().rhs, m_levels.back().sol );

  // up phase
  for( int levelIndex = numLevels - 2; levelIndex >= 0; --levelIndex )
  {
    Level const & fine = m_levels[levelIndex];
    Level const & coarse = m_levels[levelIndex + 1];
    coarse.builder->prolongation().apply( coarse.sol, fine.tmp );
    fine.sol.axpy( 1.0, fine.tmp );
    fine.matrix->residual( fine.tmp, fine.rhs, fine.rhs );
    if( fine.postsmoother )
    {
      for( integer s = 0; s < m_params.multiscale.numSmootherSweeps; ++s )
      {
        fine.postsmoother->apply( fine.rhs, fine.tmp );
        fine.sol.axpy( 1.0, fine.tmp );
        fine.matrix->residual( fine.tmp, fine.rhs, fine.rhs );
      }
    }
  }

  dst.copy( m_levels[0].sol );
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MultiscalePreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MultiscalePreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MultiscalePreconditioner< PetscInterface >;
#endif

} // namespace geosx
