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
 * @file MsrsbLevelBuilder.cpp
 */

#include "MsrsbLevelBuilder.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/MeshData.hpp"
#include "linearAlgebra/multiscale/MeshUtils.hpp"
#include "linearAlgebra/multiscale/MsrsbUtils.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
MsrsbLevelBuilder< LAI >::MsrsbLevelBuilder( string name,
                                             LinearSolverParameters::Multiscale params )
  : LevelBuilderBase< LAI >( std::move( name ), std::move( params ) ),
  m_mesh( m_name )
{}

struct BasisDescription
{
  ArrayOfSets< localIndex > supports;
  array1d< integer > supportBoundaryIndicator;
  array1d< localIndex > initialPartition;
};

/**
 * @brief Build the basic sparsity pattern for prolongation.
 *
 * Support of a coarse nodal basis function is defined as the set of fine-scale nodes
 * that are adjacent exclusively to subdomains (coarse cells or boundaries) that are
 * also adjacent to that coarse node.
 */
BasisDescription
buildNodalSupports( multiscale::MeshLevel const & fine,
                    multiscale::MeshLevel const & coarse,
                    arrayView1d< string const > const & boundaryNodeSets )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfSets< localIndex > const nodalConn =
    msrsb::buildNodalConnectivity( fine.nodeManager(), fine.cellManager() );
  ArrayOfSets< localIndex > const fineNodeToCellSubdomain =
    meshUtils::buildFineObjectToSubdomainMap( fine.nodeManager(),
                                              fine.cellManager().getExtrinsicData< meshData::CoarseCellLocalIndex >(),
                                              boundaryNodeSets );
  ArrayOfSets< localIndex > const coarseNodeToCellSubdomain =
    meshUtils::addBoundarySubdomains( coarse.nodeManager(),
                                      coarse.nodeManager().toDualRelation(),
                                      boundaryNodeSets );

  BasisDescription result;

  result.supports =
    msrsb::buildSupports( fineNodeToCellSubdomain.toViewConst(),
                          coarse.cellManager().toDualRelation().toViewConst(),
                          coarseNodeToCellSubdomain.toViewConst() );

  result.supportBoundaryIndicator =
    msrsb::findGlobalSupportBoundary( fineNodeToCellSubdomain.toViewConst() );

  result.initialPartition =
    msrsb::makeSeededPartition( nodalConn.toViewConst(),
                                coarse.nodeManager().getExtrinsicData< meshData::FineNodeLocalIndex >().toViewConst(),
                                result.supports.toViewConst() );

  return result;
}

/**
 * @brief Build the basic sparsity pattern for prolongation.
 *
 * Support of a coarse nodal basis function is defined as the set of fine-scale nodes
 * that are adjacent exclusively to subdomains (coarse cells or boundaries) that are
 * also adjacent to that coarse node.
 */
BasisDescription
buildCellSupports( multiscale::MeshLevel const & fine,
                   multiscale::MeshLevel const & coarse,
                   arrayView1d< string const > const & boundaryNodeSets )
{
  GEOSX_MARK_FUNCTION;

  // Build the nodal partition to act as dual "volumes"
  array1d< localIndex > const nodalPartition = [&]
  {
    ArrayOfSets< localIndex > const nodalConn =
      msrsb::buildNodalConnectivity( fine.nodeManager(), fine.cellManager() );
    ArrayOfSets< localIndex > const fineNodeToCellSubdomain =
      meshUtils::buildFineObjectToSubdomainMap( fine.nodeManager(),
                                                fine.cellManager().getExtrinsicData< meshData::CoarseCellLocalIndex >(),
                                                boundaryNodeSets );
    ArrayOfSets< localIndex > const coarseNodeToCellSubdomain =
      meshUtils::addBoundarySubdomains( coarse.nodeManager(),
                                        coarse.nodeManager().toDualRelation(),
                                        boundaryNodeSets );

    // Unfortunately, we need nodal supports (in order to limit nodal partition growth).
    ArrayOfSets< localIndex > const nodalSupports =
      msrsb::buildSupports( fineNodeToCellSubdomain.toViewConst(),
                            coarse.cellManager().toDualRelation(),
                            coarseNodeToCellSubdomain.toViewConst() );

    return msrsb::makeSeededPartition( nodalConn.toViewConst(),
                                       coarse.nodeManager().getExtrinsicData< meshData::FineNodeLocalIndex >(),
                                       nodalSupports.toViewConst() );
  }();

  ArrayOfSets< localIndex > const fineCellToNodalSubdomain =
    meshUtils::buildFineObjectToSubdomainMap( fine.cellManager(), nodalPartition.toViewConst(), {} );

  BasisDescription result;

  result.supports =
    msrsb::buildSupports( fineCellToNodalSubdomain.toViewConst(),
                          coarse.nodeManager().toDualRelation().toViewConst(),
                          coarse.cellManager().toDualRelation().toViewConst() );

  result.supportBoundaryIndicator =
    msrsb::findGlobalSupportBoundary( fineCellToNodalSubdomain.toViewConst() );

  result.initialPartition.resize( fine.cellManager().size() );
  result.initialPartition.setValues< parallelHostPolicy >( fine.cellManager().getExtrinsicData< meshData::CoarseCellLocalIndex >() );

  return result;
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::initializeCoarseLevel( LevelBuilderBase< LAI > & fine_level )
{
  GEOSX_MARK_FUNCTION;

  MsrsbLevelBuilder< LAI > & fine = dynamicCast< MsrsbLevelBuilder< LAI > & >( fine_level );
  m_numComp = fine.m_numComp;
  m_location = fine.m_location;

  // Coarsen the mesh
  m_mesh.buildCoarseMesh( fine.mesh(), m_params.coarsening, m_params.boundarySets );

  // Write data back to GEOSX for visualization and debug
  if( m_params.debugLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( GEOSX_FMT( "[MsRSB] {}: generated coarse grid with {} global cells and {} global nodes",
                                 m_name,
                                 m_mesh.cellManager().maxGlobalIndex() + 1,
                                 m_mesh.nodeManager().maxGlobalIndex() + 1 ) );
    GEOSX_LOG_RANK( GEOSX_FMT( "[MsRSB] {}: generated coarse grid with {} local cells and {} local nodes",
                               m_name,
                               m_mesh.cellManager().numOwnedObjects(),
                               m_mesh.nodeManager().numOwnedObjects() ) );

    m_mesh.writeCellData( { ObjectManagerBase::viewKeyStruct::ghostRankString() } );
    m_mesh.writeNodeData( { meshData::FineNodeLocalIndex::key() } );
    fine.mesh().writeCellData( { meshData::CoarseCellLocalIndex::key(),
                                 meshData::CoarseCellGlobalIndex::key() } );
    fine.mesh().writeNodeData( { meshData::CoarseNodeLocalIndex::key(),
                                 meshData::CoarseNodeGlobalIndex::key() } );
  }

  // For now, we only handle two types of basis functions - nodal and cell-centered - with specific algorithms.
  // In future this should be refactored into an extensible hierarchy of basis constructors.
  bool const isNodal = m_location == DofManager::Location::Node;
  MeshObjectManager const & coarseMgr = isNodal ? m_mesh.nodeManager() : m_mesh.cellManager();
  MeshObjectManager const & fineMgr = isNodal ? fine.mesh().nodeManager() : fine.mesh().cellManager();
  auto const buildSupports = isNodal ? buildNodalSupports : buildCellSupports;

  // Build support region definitions and tentative partition of unity
  BasisDescription const desc = buildSupports( fine.mesh(),
                                               m_mesh,
                                               m_params.boundarySets );

  // Construct global internal/boundary DoF sets
  msrsb::makeGlobalDofLists( desc.supportBoundaryIndicator,
                             m_numComp,
                             fineMgr.numOwnedObjects(),
                             fine.matrix().ilower(),
                             m_boundaryDof,
                             m_interiorDof );

  // Convert the partitioning into an actual DoF-based local matrix
  CRSMatrix< real64, globalIndex > const localProlongation =
    msrsb::buildTentativeProlongation( fineMgr,
                                       coarseMgr,
                                       desc.supports.toViewConst(),
                                       desc.initialPartition,
                                       m_numComp );

  // Assemble local pieces into a global prolongation operator and make restriction operator
  m_prolongation.create( localProlongation.toViewConst(), coarseMgr.numOwnedObjects() * m_numComp, fine.matrix().comm() );
  m_restriction = msrsb::makeRestriction( m_params, m_prolongation );

  // Create a "fake" coarse matrix (no data, just correct sizes/comms), to be computed later
  localIndex const numLocalDof = coarseMgr.numOwnedObjects() * m_numComp;
  m_matrix.createWithLocalSize( numLocalDof, numLocalDof, 0, fine.matrix().comm() );
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::initializeFineLevel( DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    string const & fieldName,
                                                    MPI_Comm const & comm )
{
  GEOSX_MARK_FUNCTION;

  m_numComp = dofManager.numComponents( fieldName );
  m_location = dofManager.location( fieldName );
  m_mesh.buildFineMesh( domain, dofManager.support( fieldName ) );

  // Create a "fake" fine matrix (no data, just correct sizes/comms for use at coarse level init)
  localIndex const numLocalDof = dofManager.numLocalDofs( fieldName );
  m_matrix.createWithLocalSize( numLocalDof, numLocalDof, 0, comm );
}

namespace
{

template< typename Matrix >
Matrix filterMatrix( Matrix const & fineMatrix,
                     integer const numComp )
{
  GEOSX_MARK_SCOPE( filter );

  // 1. Apply SC approximation
  Matrix filteredMatrix;
  fineMatrix.separateComponentFilter( filteredMatrix, numComp );

  // 2. Filter out positive off-diagonal elements (assumed positive diagonals)
  filteredMatrix.clampEntries( -LvArray::NumericLimits< real64 >::infinity, 0.0, true );

  // 3. Enforce rowsum = 0
  // 3.1. Compute rowsums
  typename Matrix::Vector rowSums;
  rowSums.create( fineMatrix.numLocalRows(), fineMatrix.comm() );
  filteredMatrix.getRowSums( rowSums, RowSumType::SumValues );

  // 3.2. Preserve Dirichlet rows by setting the diagonal update to zero
  typename Matrix::Vector diag;
  diag.create( fineMatrix.numLocalRows(), fineMatrix.comm() );
  filteredMatrix.extractDiagonal( diag );
  forAll< parallelHostPolicy >( diag.localSize(), [diagData = diag.values(),
                                                   rowSumData = rowSums.open()]( localIndex const localRow )
  {
    if( isEqual( diagData[localRow], rowSumData[localRow] ) )
    {
      rowSumData[localRow] = 0.0;
    }
  } );
  rowSums.close();

  // 3.3. Subtract the nonzero rowsums from diagonal elements
  filteredMatrix.addDiagonal( rowSums, -1.0 );

  return filteredMatrix;
}

template< typename MATRIX >
auto makeJacobiMatrix( MATRIX && fineMatrix,
                       real64 const omega )
{
  GEOSX_MARK_SCOPE( jacobi );
  using Matrix = std::remove_const_t< TYPEOFREF( fineMatrix ) >;

  // 0. Copy or move input matrix into a new object
  Matrix iterMatrix( std::forward< MATRIX >( fineMatrix ) );

  // 1. Compute -w * Dinv * A;
  typename Matrix::Vector diag;
  diag.create( iterMatrix.numLocalRows(), iterMatrix.comm() );
  iterMatrix.extractDiagonal( diag );
  diag.reciprocal();
  diag.scale( -omega );
  iterMatrix.leftScale( diag );

  // 2. Compute I - w * Dinv * A by adding identity diagonal
  diag.set( 1.0 );
  iterMatrix.addDiagonal( diag, 1.0 );
  return iterMatrix;
}

template< typename Matrix >
Matrix makeIterationMatrix( Matrix const & fineMatrix,
                            integer const numComp,
                            real64 const omega,
                            integer const debugLevel,
                            string const & debugPrefix )
{
  GEOSX_MARK_SCOPE( make_iter_matrix );

  Matrix filteredMatrix = filterMatrix( fineMatrix, numComp );
  if( debugLevel >= 4 )
  {
    filteredMatrix.write( debugPrefix + "_filtered.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  Matrix jacobiMatrix = makeJacobiMatrix( std::move( filteredMatrix ), omega );
  if( debugLevel >= 4 )
  {
    jacobiMatrix.write( debugPrefix + "_jacobi.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  return jacobiMatrix;
}

template< typename Matrix >
integer iterateBasis( Matrix const & jacobiMatrix,
                      arrayView1d< globalIndex const > const & boundaryDof,
                      arrayView1d< globalIndex const > const & interiorDof,
                      integer const maxIter,
                      real64 const tolerance,
                      integer const checkFreq,
                      integer const debugLevel,
                      string const & name,
                      Matrix & prolongation )
{
  GEOSX_MARK_SCOPE( iterate );

  auto const saveForDebug = [&]( string const & suffix, integer const minDebugLevel )
  {
    if( debugLevel >= minDebugLevel )
    {
      GEOSX_MARK_SCOPE( writeProlongationMatrix );
      prolongation.write( GEOSX_FMT( "{}_P_{}.mtx", name, suffix ), LAIOutputFormat::MATRIX_MARKET );
    }
  };

  Matrix P( prolongation );
  integer iter = 0;
  real64 norm = LvArray::NumericLimits< real64 >::max;

  auto const computeAndLogConvergenceNorm = [&]()
  {
    GEOSX_MARK_SCOPE( check );
    P.addEntries( prolongation, MatrixPatternOp::Same, -1.0 );
    norm = P.normMax( interiorDof );
    GEOSX_LOG_RANK_0_IF( debugLevel >= 3, GEOSX_FMT( "[MsRSB] {}: iter = {}, conv = {:e}", name, iter, norm ) );
  };

  saveForDebug( "init", 4 );
  while( iter < maxIter && norm > tolerance )
  {
    // Keep 1-based iteration index for convenience
    ++iter;

    // Perform a step of Jacobi
    Matrix Ptemp;
    {
      GEOSX_MARK_SCOPE( multiply );
      jacobiMatrix.multiply( prolongation, Ptemp );
    }

    // Restrict to the predefined prolongation pattern
    {
      GEOSX_MARK_SCOPE( restrict );
      P.zero();
      P.addEntries( Ptemp, MatrixPatternOp::Restrict, 1.0 );
    }

    // Rescale to preserve partition of unity
    {
      GEOSX_MARK_SCOPE( rescale );
      P.rescaleRows( boundaryDof, RowSumType::SumValues );
    }

    // Switch over to new prolongation operator
    std::swap( P, prolongation );
    saveForDebug( std::to_string( iter ), 6 );

    // Compute update norm, check convergence
    if( iter % checkFreq == 0 )
    {
      computeAndLogConvergenceNorm();
    }
  }

  // Compute update norm and check convergence one final time if needed (in case we ran out of iterations)
  if( iter % checkFreq != 0 )
  {
    computeAndLogConvergenceNorm();
  }

  GEOSX_LOG_RANK_0_IF( debugLevel >= 1, GEOSX_FMT( "[MsRSB] {}: {} in {} iterations", name, norm <= tolerance ? "converged" : "failed to converge", iter ) );

  saveForDebug( "conv", 4 );
  return iter;
}

} // namespace

template< typename LAI >
void MsrsbLevelBuilder< LAI >::compute( Matrix const & fineMatrix )
{
  GEOSX_MARK_FUNCTION;

  // Compute prolongation
  GEOSX_LOG_RANK_0_IF( m_params.debugLevel >= 2, GEOSX_FMT( "[MsRSB] {}: building iteration matrix", m_name ) );
  Matrix const jacobiMatrix = makeIterationMatrix( fineMatrix,
                                                   m_numComp,
                                                   m_params.msrsb.relaxation,
                                                   m_params.debugLevel,
                                                   m_name );

  GEOSX_LOG_RANK_0_IF( m_params.debugLevel >= 2, GEOSX_FMT( "[MsRSB] {}: performing basis iteration", m_name ) );
  m_lastNumIter = iterateBasis( jacobiMatrix,
                                m_boundaryDof,
                                m_interiorDof,
                                m_params.msrsb.maxIter,
                                m_params.msrsb.tolerance,
                                m_lastNumIter <= m_params.msrsb.checkFrequency ? 1 : m_params.msrsb.checkFrequency,
                                m_params.debugLevel,
                                m_name,
                                m_prolongation );

  if( m_lastNumIter > 1 && m_params.debugLevel >= 5 )
  {
    GEOSX_MARK_SCOPE( plotBasis );
    writeProlongationForDebug();
  }

  // Recompute coarse operator - only if prolongation took a nontrivial number of iterations to converge
  if( m_lastNumIter > 1 || !m_matrix.ready() )
  {
    GEOSX_MARK_SCOPE( RAP );
    GEOSX_LOG_RANK_0_IF( m_params.debugLevel >= 2, GEOSX_FMT( "[MsRSB] {}: computing RAP", m_name ) );
    if( m_params.galerkin )
    {
      fineMatrix.multiplyPtAP( m_prolongation, m_matrix );
    }
    else
    {
      Matrix const & restriction = dynamicCast< Matrix const & >( *m_restriction );
      fineMatrix.multiplyRAP( m_prolongation, restriction, m_matrix );
    }
  }

  if( m_params.debugLevel >= 4 )
  {
    GEOSX_MARK_SCOPE( writeLevelMatrix );
    m_matrix.write( m_name + ".mtx", LAIOutputFormat::MATRIX_MARKET );
  }
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::writeProlongationForDebug() const
{
  if( m_location == DofManager::Location::Node )
  {
    msrsb::writeProlongation( m_prolongation,
                              m_name,
                              m_numComp,
                              *m_mesh.fineMesh()->domain(),
                              m_mesh.fineMesh()->nodeManager(),
                              [&]( std::vector< string > const & names )
                              { m_mesh.fineMesh()->writeNodeData( names ); } );
  }
  else
  {
    msrsb::writeProlongation( m_prolongation,
                              m_name,
                              m_numComp,
                              *m_mesh.fineMesh()->domain(),
                              m_mesh.fineMesh()->cellManager(),
                              [&]( std::vector< string > const & names )
    { m_mesh.fineMesh()->writeCellData( names ); } );
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MsrsbLevelBuilder< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsrsbLevelBuilder< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsrsbLevelBuilder< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geosx
