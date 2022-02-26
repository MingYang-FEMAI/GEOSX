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
 * @file MsrsbLevelBuilder.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP

#include "linearAlgebra/multiscale/MeshLevel.hpp"
#include "linearAlgebra/multiscale/LevelBuilderBase.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
class MsrsbLevelBuilder : public LevelBuilderBase< LAI >
{
public:

  /// Alias for base type
  using Base = LevelBuilderBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Alias for operator type
  using Operator = typename Base::Operator;

  explicit MsrsbLevelBuilder( string name, LinearSolverParameters::Multiscale params );

  virtual Operator const & prolongation() const override
  {
    return m_prolongation;
  }

  virtual Operator const & restriction() const override
  {
    return *m_restriction;
  }

  virtual Matrix const & matrix() const override
  {
    return m_matrix;
  }

  virtual void initializeFineLevel( DomainPartition & domain,
                                    DofManager const & dofManager,
                                    string const & fieldName,
                                    MPI_Comm const & comm ) override;

  virtual void initializeCoarseLevel( LevelBuilderBase< LAI > & fine_level ) override;

  virtual void compute( Matrix const & fineMatrix ) override;

  multiscale::MeshLevel & mesh() { return m_mesh; }
  multiscale::MeshLevel const & mesh() const { return m_mesh; }

private:

  void writeProlongationForDebug() const;

  using Base::m_params;
  using Base::m_name;

  /// Number of dof components
  integer m_numComp = 1;

  /// Dof location (cell or node)
  DofManager::Location m_location = DofManager::Location::Node;

  /// Mesh description at current level
  multiscale::MeshLevel m_mesh;

  /// Prolongation matrix P
  Matrix m_prolongation;

  /// Restriction (kept as abstract operator to allow for memory efficiency, e.g. when R = P^T)
  std::unique_ptr< Operator > m_restriction;

  /// Level operator matrix
  Matrix m_matrix;

  /// List of nodes on global boundary
  array1d< globalIndex > m_boundaryDof;

  /// List of nodes that are interior
  array1d< globalIndex > m_interiorDof;

  /// Previous number of smoothing iterations
  integer m_lastNumIter = std::numeric_limits< integer >::max();
};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP
