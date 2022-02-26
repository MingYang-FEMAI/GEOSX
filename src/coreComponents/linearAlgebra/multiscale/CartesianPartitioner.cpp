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
 * @file CartesianPartitioner.cpp
 */

#include "CartesianPartitioner.hpp"

#include "linearAlgebra/multiscale/MeshData.hpp"

namespace geosx
{
namespace multiscale
{

CartesianPartitioner::CartesianPartitioner( LinearSolverParameters::Multiscale::Coarsening params )
  : m_params( std::move( params ) )
{}

localIndex CartesianPartitioner::generate( MeshLevel const & mesh,
                                           arrayView1d< localIndex > const & partition )
{
  GEOSX_MARK_FUNCTION;

  MeshObjectManager const & cellManager = mesh.cellManager();

  GEOSX_THROW_IF_NE_MSG( m_params.ratio.size(), 3,
                         "Cartesian partitioning requires 3 coarsening ratio values",
                         InputError );
  GEOSX_THROW_IF( !cellManager.hasExtrinsicData< meshData::CartesianIndex >(),
                  "Cartesian partitioning is only compatible with InternalMesh",
                  InputError );

  // Trivial case when the rank doesn't have any cells in the domain of multiscale
  if( cellManager.numOwnedObjects() == 0 )
  {
    return 0;
  }

  arrayView2d< integer const > const cartIndex = cellManager.getExtrinsicData< meshData::CartesianIndex >().toViewConst();

  // Find the range of locally owned cell cartesian indices
  integer loCartIndex[3]{};
  integer hiCartIndex[3]{};
  for( int dim = 0; dim < 3; ++dim )
  {
    RAJA::ReduceMin< parallelHostReduce, integer > lo( std::numeric_limits< integer >::max() );
    RAJA::ReduceMax< parallelHostReduce, integer > hi( std::numeric_limits< integer >::min() );
    forAll< parallelHostPolicy >( cellManager.numOwnedObjects(), [=]( localIndex const i )
    {
      if( cartIndex[i][dim] >= 0 )
      {
        lo.min( cartIndex[i][dim] );
        hi.max( cartIndex[i][dim] );
      }
    } );
    loCartIndex[dim] = lo;
    hiCartIndex[dim] = hi;
  }

  // Sanity check
  GEOSX_ASSERT_GE( hiCartIndex[0], 0 );
  GEOSX_ASSERT_GE( hiCartIndex[1], 0 );
  GEOSX_ASSERT_GE( hiCartIndex[2], 0 );

  // Compute cartesian sizes of coarse grid
  integer cartRatio[3];
  for( int dim = 0; dim < 3; ++dim )
  {
    integer const cartSizeFine = hiCartIndex[dim] - loCartIndex[dim] + 1;
    integer const cartSizeCoarse = integer( std::ceil( cartSizeFine / m_params.ratio[dim] ) ); // rounded up
    m_cartNumPart[dim] = cartSizeCoarse;
    cartRatio[dim] = ( cartSizeFine + cartSizeCoarse - 1) / cartSizeCoarse; // rounded up
  }

  localIndex const cartPartStride[3] = { 1, m_cartNumPart[0], m_cartNumPart[0] * m_cartNumPart[1] };

  // Compute cartesian coarse cell indices
  forAll< parallelHostPolicy >( cellManager.numOwnedObjects(), [=]( localIndex const i )
  {
    localIndex part = 0;
    for( int dim = 0; dim < 3; ++dim )
    {
      integer const cartIndexCoarse = ( cartIndex[i][dim] - loCartIndex[dim] ) / cartRatio[dim];
      part += cartIndexCoarse * cartPartStride[dim];
    }
    partition[i] = part;
  } );

  return m_cartNumPart[0] * m_cartNumPart[1] * m_cartNumPart[2];
}

void CartesianPartitioner::setCoarseData( multiscale::MeshLevel & coarseMesh ) const
{
  array2d< integer > & cartIndex =
    coarseMesh.cellManager().registerExtrinsicData< meshData::CartesianIndex >( {} ).reference();
  cartIndex.resizeDimension< 1 >( 3 );

  localIndex const cartPartStride[3] = { 1, m_cartNumPart[0], m_cartNumPart[0] * m_cartNumPart[1] };

  forAll< parallelHostPolicy >( coarseMesh.cellManager().numOwnedObjects(),
                                [cartIndex = cartIndex.toView(), cartPartStride]( localIndex const i )
  {
    localIndex cellIndex = i;
    for( int dim = 2; dim >= 0; --dim )
    {
      cartIndex[i][dim] = LvArray::integerConversion< integer >( cellIndex / cartPartStride[dim] );
      cellIndex = cellIndex % cartPartStride[dim];
    }
  } );
}

} // namespace multiscale
} // namespace geosx
