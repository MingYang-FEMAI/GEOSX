/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverParameters.cpp
 */

#include "LinearSolverParameters.hpp"

namespace geosx
{

using namespace dataRepository;

class MetisParametersInput : public dataRepository::Group
{
public:

  /// Constructor
  MetisParametersInput( string const & name,
                        Group * const parent,
                        LinearSolverParameters::Multiscale::Coarsening::Metis & params );

  virtual Group * createChild( string const & childKey, string const & childName ) override final
  {
    GEOSX_UNUSED_VAR( childKey, childName );
    return nullptr;
  }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr char const * methodString()    { return "method"; }
    static constexpr char const * minCommonNodesString() { return "minCommonNodes"; }
    static constexpr char const * ufactorString()   { return "ufactor"; }
  };

private:

  LinearSolverParameters::Multiscale::Coarsening::Metis & m_parameters;
};

MetisParametersInput::MetisParametersInput( string const & name,
                                            Group * const parent,
                                            LinearSolverParameters::Multiscale::Coarsening::Metis & params )
  :
  Group( name, parent ),
  m_parameters( params )
{
  registerWrapper( viewKeyStruct::methodString(), &m_parameters.method ).
    setApplyDefaultValue( m_parameters.method ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "METIS partitioning method, one of: "
                    "``" + EnumStrings< LinearSolverParameters::Multiscale::Coarsening::Metis::Method >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::minCommonNodesString(), &m_parameters.minCommonNodes ).
    setApplyDefaultValue( m_parameters.minCommonNodes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum number of nodes shared between two cells when constructing the connectivity graph" );

  registerWrapper( viewKeyStruct::ufactorString(), &m_parameters.ufactor ).
    setApplyDefaultValue( m_parameters.ufactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "METIS ufactor parameter, affects partitioning balance/edgecut tradeoff" );
}

class CoarseningParametersInput : public dataRepository::Group
{
public:

  /// Constructor
  CoarseningParametersInput( string const & name,
                             Group * const parent,
                             LinearSolverParameters::Multiscale::Coarsening & params );

  virtual Group * createChild( string const & childKey, string const & childName ) override final
  {
    GEOSX_UNUSED_VAR( childKey, childName );
    return nullptr;
  }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr char const * partitionTypeString() { return "partitionType"; }
    static constexpr char const * ratioString()         { return "ratio"; }
    static constexpr char const * minLocalDofString()   { return "minLocalDof"; }
    static constexpr char const * minGlobalDofString()  { return "minGlobalDof"; }
  };

  /// Keys appearing in XML
  struct groupKeyStruct
  {
    static constexpr char const * metisString() { return "Metis"; }
  };

private:

  LinearSolverParameters::Multiscale::Coarsening & m_parameters;
};

CoarseningParametersInput::CoarseningParametersInput( string const & name,
                                                      Group * const parent,
                                                      LinearSolverParameters::Multiscale::Coarsening & params )
  :
  Group( name, parent ),
  m_parameters( params )
{
  registerWrapper( viewKeyStruct::partitionTypeString(), &m_parameters.partitionType ).
    setApplyDefaultValue( m_parameters.partitionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Partition type for generating coarse aggregates. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Multiscale::Coarsening::PartitionType >::concat( "``, ``" ) + "``" );

  array1d< real64 > & coarseningRatio = registerWrapper( viewKeyStruct::ratioString(), &m_parameters.ratio ).
                                          setInputFlag( InputFlags::OPTIONAL ).
                                          setDescription( "Coarsening ratio (number of fine cells per coarse cell)" ).reference();
  coarseningRatio.resize( 3 );
  coarseningRatio[0] = 8;
  coarseningRatio[1] = 8;
  coarseningRatio[2] = 8;

  registerWrapper( viewKeyStruct::minLocalDofString(), &m_parameters.minLocalDof ).
    setApplyDefaultValue( m_parameters.minLocalDof ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Limit of coarsening on current rank (i.e. keep a local coarsening ratio of 1 once this problem size reached)" );

  registerWrapper( viewKeyStruct::minGlobalDofString(), &m_parameters.minGlobalDof ).
    setApplyDefaultValue( m_parameters.minGlobalDof ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "limit of coarsening across all ranks (i.e. trim the grid hierarchy globally)" );

  registerGroup( groupKeyStruct::metisString(),
                 std::make_unique< MetisParametersInput >( groupKeyStruct::metisString(), this, m_parameters.metis ) ).
    setInputFlags( InputFlags::OPTIONAL );
}


class MsrsbParametersInput : public dataRepository::Group
{
public:

  /// Constructor
  MsrsbParametersInput( string const & name,
                        Group * const parent,
                        LinearSolverParameters::Multiscale::MsRSB & params );

  virtual Group * createChild( string const & childKey, string const & childName ) override final
  {
    GEOSX_UNUSED_VAR( childKey, childName );
    return nullptr;
  }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr char const * maxIterString()            { return "maxIter"; }
    static constexpr char const * toleranceString()          { return "tolerance"; }
    static constexpr char const * relaxationString()         { return "relaxation"; }
    static constexpr char const * checkFrequencyString()     { return "checkFrequency"; }
  };

private:

  LinearSolverParameters::Multiscale::MsRSB & m_parameters;
};

MsrsbParametersInput::MsrsbParametersInput( string const & name,
                                            Group * const parent,
                                            LinearSolverParameters::Multiscale::MsRSB & params )
  :
  Group( name, parent ),
  m_parameters( params )
{
  registerWrapper( viewKeyStruct::maxIterString(), &m_parameters.maxIter ).
    setApplyDefaultValue( m_parameters.maxIter ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of MsRSB basis smoothing iterations" );

  registerWrapper( viewKeyStruct::toleranceString(), &m_parameters.tolerance ).
    setApplyDefaultValue( m_parameters.tolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "MsRSB basis smoothing iteration tolerance" );

  registerWrapper( viewKeyStruct::relaxationString(), &m_parameters.relaxation ).
    setApplyDefaultValue( m_parameters.relaxation ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "MsRSB basis smoothing iteration relaxation parameter" );

  registerWrapper( viewKeyStruct::checkFrequencyString(), &m_parameters.checkFrequency ).
    setApplyDefaultValue( m_parameters.checkFrequency ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "MsRSB basis smoothing convergence check frequency" );
}


class MultiscaleParametersInput : public dataRepository::Group
{
public:

  /// Constructor
  MultiscaleParametersInput( string const & name,
                             Group * const parent,
                             LinearSolverParameters::Multiscale & params );

  virtual Group * createChild( string const & childKey, string const & childName ) override final
  {
    GEOSX_UNUSED_VAR( childKey, childName );
    return nullptr;
  }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr char const * basisTypeString()               { return "basisType"; }
    static constexpr char const * maxLevelsString()               { return "maxLevels"; }
    static constexpr char const * numSmootherSweepsString()       { return "numSmootherSweeps"; }
    static constexpr char const * preOrPostSmoothingString()      { return "preOrPostSmoothing"; }
    static constexpr char const * smootherTypeString()            { return "smootherType"; }
    static constexpr char const * boundarySets()                  { return "boundarySets"; }
    static constexpr char const * debugLevel()                    { return "debugLevel"; }
    static constexpr char const * coarseTypeString()              { return "coarseType"; }
  };

  /// Keys appearing in XML
  struct groupKeyStruct
  {
    static constexpr char const * coarseningString() { return "Coarsening"; }
    static constexpr char const * msrsbString() { return "MsRSB"; }
  };

private:

  LinearSolverParameters::Multiscale & m_parameters;
};

MultiscaleParametersInput::MultiscaleParametersInput( string const & name,
                                                      Group * const parent,
                                                      LinearSolverParameters::Multiscale & params )
  :
  Group( name, parent ),
  m_parameters( params )
{
  registerWrapper( viewKeyStruct::basisTypeString(), &m_parameters.basisType ).
    setApplyDefaultValue( m_parameters.basisType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Multiscale basis type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Multiscale::BasisType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::maxLevelsString(), &m_parameters.maxLevels ).
    setApplyDefaultValue( m_parameters.maxLevels ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of multiscale grid levels (including fine)" );

  registerWrapper( viewKeyStruct::numSmootherSweepsString(), &m_parameters.numSmootherSweeps ).
    setApplyDefaultValue( m_parameters.numSmootherSweeps ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of smoother sweeps" );

  registerWrapper( viewKeyStruct::preOrPostSmoothingString(), &m_parameters.preOrPostSmoothing ).
    setApplyDefaultValue( m_parameters.preOrPostSmoothing ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Pre and/or post smoothing (``" + EnumStrings< LinearSolverParameters::AMG::PreOrPost >::concat( "``, ``" ) + "``)" );

  registerWrapper( viewKeyStruct::smootherTypeString(), &m_parameters.smootherType ).
    setApplyDefaultValue( m_parameters.smootherType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Smoother type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::PreconditionerType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::boundarySets(), &m_parameters.boundarySets ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of node set names that denote global domain boundaries, improves interpolation when provided." );

  registerWrapper( viewKeyStruct::debugLevel(), &m_parameters.debugLevel ).
    setApplyDefaultValue( m_parameters.debugLevel ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Debug level (0 - no debug, 1 - basic progress messages, 2 - detailed output and matrix dumps)" );

  registerWrapper( viewKeyStruct::coarseTypeString(), &m_parameters.coarseType ).
    setApplyDefaultValue( m_parameters.coarseType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coarsest level solver type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::PreconditionerType >::concat( "``, ``" ) + "``" );

  registerGroup( groupKeyStruct::coarseningString(),
                 std::make_unique< CoarseningParametersInput >( groupKeyStruct::coarseningString(), this, m_parameters.coarsening ) ).
    setInputFlags( InputFlags::OPTIONAL );

  registerGroup( groupKeyStruct::msrsbString(),
                 std::make_unique< MsrsbParametersInput >( groupKeyStruct::msrsbString(), this, m_parameters.msrsb ) ).
    setInputFlags( InputFlags::OPTIONAL );
}

LinearSolverParametersInput::LinearSolverParametersInput( string const & name,
                                                          Group * const parent )
  :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solverTypeString(), &m_parameters.solverType ).
    setApplyDefaultValue( m_parameters.solverType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Linear solver type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::SolverType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::preconditionerTypeString(), &m_parameters.preconditionerType ).
    setApplyDefaultValue( m_parameters.preconditionerType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Preconditioner type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::PreconditionerType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::stopIfErrorString(), &m_parameters.stopIfError ).
    setApplyDefaultValue( m_parameters.stopIfError ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to stop the simulation if the linear solver reports an error" );

  registerWrapper( viewKeyStruct::directCheckResidualString(), &m_parameters.direct.checkResidual ).
    setApplyDefaultValue( m_parameters.direct.checkResidual ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to check the linear system solution residual" );

  registerWrapper( viewKeyStruct::directEquilString(), &m_parameters.direct.equilibrate ).
    setApplyDefaultValue( m_parameters.direct.equilibrate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to scale the rows and columns of the matrix" );

  registerWrapper( viewKeyStruct::directColPermString(), &m_parameters.direct.colPerm ).
    setApplyDefaultValue( m_parameters.direct.colPerm ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How to permute the columns. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Direct::ColPerm >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::directRowPermString(), &m_parameters.direct.rowPerm ).
    setApplyDefaultValue( m_parameters.direct.rowPerm ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How to permute the rows. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Direct::RowPerm >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::directReplTinyPivotString(), &m_parameters.direct.replaceTinyPivot ).
    setApplyDefaultValue( m_parameters.direct.replaceTinyPivot ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to replace tiny pivots by sqrt(epsilon)*norm(A)" );

  registerWrapper( viewKeyStruct::directIterRefString(), &m_parameters.direct.iterativeRefine ).
    setApplyDefaultValue( m_parameters.direct.iterativeRefine ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to perform iterative refinement" );

  registerWrapper( viewKeyStruct::directParallelString(), &m_parameters.direct.parallel ).
    setApplyDefaultValue( m_parameters.direct.parallel ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to use a parallel solver (instead of a serial one)" );

  registerWrapper( viewKeyStruct::krylovMaxIterString(), &m_parameters.krylov.maxIterations ).
    setApplyDefaultValue( m_parameters.krylov.maxIterations ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum iterations allowed for an iterative solver" );

  registerWrapper( viewKeyStruct::krylovMaxRestartString(), &m_parameters.krylov.maxRestart ).
    setApplyDefaultValue( m_parameters.krylov.maxRestart ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum iterations before restart (GMRES only)" );

  registerWrapper( viewKeyStruct::krylovTolString(), &m_parameters.krylov.relTolerance ).
    setApplyDefaultValue( m_parameters.krylov.relTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Relative convergence tolerance of the iterative method\n"
                    "If the method converges, the iterative solution :math:`\\mathsf{x}_k` is such that\n"
                    "the relative residual norm satisfies:\n"
                    ":math:`\\left\\lVert \\mathsf{b} - \\mathsf{A} \\mathsf{x}_k \\right\\rVert_2` < ``" +
                    string( viewKeyStruct::krylovTolString() ) + "`` * :math:`\\left\\lVert\\mathsf{b}\\right\\rVert_2`" );

  registerWrapper( viewKeyStruct::krylovAdaptiveTolString(), &m_parameters.krylov.useAdaptiveTol ).
    setApplyDefaultValue( m_parameters.krylov.useAdaptiveTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use Eisenstat-Walker adaptive linear tolerance" );

  registerWrapper( viewKeyStruct::krylovWeakTolString(), &m_parameters.krylov.weakestTol ).
    setApplyDefaultValue( m_parameters.krylov.weakestTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Weakest-allowed tolerance for adaptive method" );

  registerWrapper( viewKeyStruct::amgNumSweepsString(), &m_parameters.amg.numSweeps ).
    setApplyDefaultValue( m_parameters.amg.numSweeps ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG smoother sweeps" );

  registerWrapper( viewKeyStruct::amgSmootherString(), &m_parameters.amg.smootherType ).
    setApplyDefaultValue( m_parameters.amg.smootherType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG smoother type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::SmootherType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::amgCoarseString(), &m_parameters.amg.coarseType ).
    setApplyDefaultValue( m_parameters.amg.coarseType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG coarsest level solver/smoother type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::CoarseType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::amgCoarseningString(), &m_parameters.amg.coarseningType ).
    setApplyDefaultValue( "HMIS" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG coarsening algorithm\n"
                    "Available options are: TODO" );

  registerWrapper( viewKeyStruct::amgInterpolationString(), &m_parameters.amg.interpolationType ).
    setApplyDefaultValue( 6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG interpolation algorithm\n"
                    "Available options are: TODO" );

  registerWrapper( viewKeyStruct::amgNumFunctionsString(), &m_parameters.amg.numFunctions ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG number of functions\n"
                    "Available options are: TODO" );

  registerWrapper( viewKeyStruct::amgAggresiveNumLevelsString(), &m_parameters.amg.aggresiveNumLevels ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG number levels for aggressive coarsening \n"
                    "Available options are: TODO" );

  registerWrapper( viewKeyStruct::amgThresholdString(), &m_parameters.amg.threshold ).
    setApplyDefaultValue( m_parameters.amg.threshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG strength-of-connection threshold" );

  registerWrapper( viewKeyStruct::amgNullSpaceTypeString(), &m_parameters.amg.nullSpaceType ).
    setApplyDefaultValue( m_parameters.amg.nullSpaceType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG near null space approximation. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::NullSpaceType >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeyStruct::iluFillString(), &m_parameters.ifact.fill ).
    setApplyDefaultValue( m_parameters.ifact.fill ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "ILU(K) fill factor" );

  registerWrapper( viewKeyStruct::iluThresholdString(), &m_parameters.ifact.threshold ).
    setApplyDefaultValue( m_parameters.ifact.threshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "ILU(T) threshold factor" );

  registerGroup( groupKeyStruct::multiscaleString(),
                 std::make_unique< MultiscaleParametersInput >( groupKeyStruct::multiscaleString(), this, m_parameters.multiscale ) ).
    setInputFlags( InputFlags::OPTIONAL );
}

void LinearSolverParametersInput::postProcessInput()
{
  m_parameters.logLevel = getLogLevel();

  static const std::set< integer > binaryOptions = { 0, 1 };

  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.stopIfError ) == 0, viewKeyStruct::stopIfErrorString() << " option can be either 0 (false) or 1 (true)" );
  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.direct.checkResidual ) == 0, viewKeyStruct::directCheckResidualString() << " option can be either 0 (false) or 1 (true)" );
  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.direct.equilibrate ) == 0, viewKeyStruct::directEquilString() << " option can be either 0 (false) or 1 (true)" );
  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.direct.replaceTinyPivot ) == 0, viewKeyStruct::directReplTinyPivotString() << " option can be either 0 (false) or 1 (true)" );
  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.direct.iterativeRefine ) == 0, viewKeyStruct::directIterRefString() << " option can be either 0 (false) or 1 (true)" );
  GEOSX_ERROR_IF( binaryOptions.count( m_parameters.direct.parallel ) == 0, viewKeyStruct::directParallelString() << " option can be either 0 (false) or 1 (true)" );

  GEOSX_ERROR_IF_LT_MSG( m_parameters.krylov.maxIterations, 0, "Invalid value of " << viewKeyStruct::krylovMaxIterString() );
  GEOSX_ERROR_IF_LT_MSG( m_parameters.krylov.maxRestart, 0, "Invalid value of " << viewKeyStruct::krylovMaxRestartString() );

  GEOSX_ERROR_IF_LT_MSG( m_parameters.krylov.relTolerance, 0.0, "Invalid value of " << viewKeyStruct::krylovTolString() );
  GEOSX_ERROR_IF_GT_MSG( m_parameters.krylov.relTolerance, 1.0, "Invalid value of " << viewKeyStruct::krylovTolString() );

  GEOSX_ERROR_IF_LT_MSG( m_parameters.ifact.fill, 0, "Invalid value of " << viewKeyStruct::iluFillString() );
  GEOSX_ERROR_IF_LT_MSG( m_parameters.ifact.threshold, 0.0, "Invalid value of " << viewKeyStruct::iluThresholdString() );

  GEOSX_ERROR_IF_LT_MSG( m_parameters.amg.numSweeps, 0, "Invalid value of " << viewKeyStruct::amgNumSweepsString() );
  GEOSX_ERROR_IF_LT_MSG( m_parameters.amg.threshold, 0.0, "Invalid value of " << viewKeyStruct::amgThresholdString() );
  GEOSX_ERROR_IF_GT_MSG( m_parameters.amg.threshold, 1.0, "Invalid value of " << viewKeyStruct::amgThresholdString() );

  // TODO input validation for other AMG parameters ?
}

Group * LinearSolverParametersInput::createChild( string const & childKey,
                                                  string const & childName )
{
  GEOSX_UNUSED_VAR( childKey, childName );
  return nullptr;
}

REGISTER_CATALOG_ENTRY( Group, LinearSolverParametersInput, string const &, Group * const )

} // namespace geosx
