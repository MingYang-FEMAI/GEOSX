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
 * @file MsrsbUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/multiscale/MeshObjectManager.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{
namespace multiscale
{

class MeshObjectManager;

namespace msrsb
{

ArrayOfSets< localIndex >
buildNodalConnectivity( MeshObjectManager const & nodeManager,
                        MeshObjectManager const & cellManager );

array1d< localIndex >
makeSeededPartition( ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & seeds,
                     ArrayOfSetsView< localIndex const > const & supports );

/**
 * @brief Build a map of support regions through subdomain adjacency.
 * @param fineObjectToSubdomain map of fine-scale objects (points) to adjacent subdomains
 *        ("virtual" boundary subdomains, if included, must have negative indices)
 * @param subdomainToCoarseObject map of subdomains to adjacent coarse-scale objects (points).
 *        Does not need to include boundary subdomains.
 * @param coarseObjectToSubdomain map of coarse-scale objects (points) to adjacent subdomains
 *        (generally a transpose of @p subdomainToCoarseObject)
 * @param supportBoundaryIndicator auxiliary output array that marks points on global support boundary with a nonzero value (1).
 * @return a map of fine objects (points) to a list of coarse objects (points) to the support of which they belong
 */
ArrayOfSets< localIndex >
buildSupports( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain,
               ArrayOfSetsView< localIndex const > const & subdomainToCoarseObject,
               ArrayOfSetsView< localIndex const > const & coarseObjectToSubdomain );

array1d< integer >
findGlobalSupportBoundary( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain );

SparsityPattern< globalIndex >
buildProlongationSparsity( MeshObjectManager const & fineManager,
                           MeshObjectManager const & coarseManager,
                           ArrayOfSetsView< localIndex const > const & supports,
                           integer const numComp );

CRSMatrix< real64, globalIndex >
buildTentativeProlongation( MeshObjectManager const & fineManager,
                            MeshObjectManager const & coarseManager,
                            ArrayOfSetsView< localIndex const > const & supports,
                            arrayView1d< localIndex const > const & initPart,
                            integer const numComp );

void makeGlobalDofLists( arrayView1d< integer const > const & indicator,
                         integer const numComp,
                         localIndex const numLocalObjects,
                         globalIndex const firstLocalDof,
                         array1d< globalIndex > & boundaryDof,
                         array1d< globalIndex > & interiorDof );

template< typename MATRIX >
std::unique_ptr< LinearOperator< typename MATRIX::Vector > >
makeRestriction( LinearSolverParameters::Multiscale const & params,
                 MATRIX const & prolongation )
{
  std::unique_ptr< LinearOperator< typename MATRIX::Vector > > restriction;
  if( params.galerkin )
  {
    // Make a transpose operator with a reference to P, which will be computed later
    restriction = std::make_unique< TransposeOperator< MATRIX > >( prolongation );
  }
  else
  {
    // Make an explicit transpose of tentative (initial) P
    MATRIX R;
    prolongation.transpose( R );
    restriction = std::make_unique< MATRIX >( std::move( R ) );
  }
  return restriction;
}

template< typename Matrix >
void writeProlongation( Matrix const & prolongation,
                        string const & prefix,
                        integer const numComp,
                        DomainPartition & domain,
                        multiscale::MeshObjectManager & fineManager,
                        std::function< void ( std::vector< string > const & ) > const & writeFunc )
{
  std::vector< string > bNames{ "X ", "Y ", "Z " };
  std::vector< string > cNames{ " x", " y", " z" };

  globalIndex const numCoarsePoints = prolongation.numGlobalCols() / numComp;
  int const labelWidth = static_cast< int >( std::log10( numCoarsePoints ) ) + 1;

  std::vector< arrayView3d< real64 > > views;
  std::vector< string > names;

  for( globalIndex cidx = 0; cidx < numCoarsePoints; ++cidx )
  {
    string const name = GEOSX_FMT( "{}_P_{:{}}", prefix, cidx, labelWidth );
    auto & wrapper = fineManager.registerWrapper< array3d< real64 > >( name ).
                       setDimLabels( 1, { bNames.begin(), bNames.begin() + numComp } ).
                       setDimLabels( 2, { cNames.begin(), cNames.begin() + numComp } ).
                       setPlotLevel( dataRepository::PlotLevel::LEVEL_0 );
    wrapper.reference().resize( fineManager.size(), numComp, numComp );
    views.push_back( wrapper.referenceAsView() );
    names.push_back( name );
  }

  array1d< globalIndex > colIndices;
  array1d< real64 > values;

  for( localIndex localRow = 0; localRow < prolongation.numLocalRows(); ++localRow )
  {
    globalIndex const globalRow = prolongation.ilower() + localRow;
    localIndex const rowObj = localRow / numComp;
    integer const rowComp = static_cast< integer >( localRow % numComp );

    localIndex const numValues = prolongation.rowLength( globalRow );
    colIndices.resize( numValues );
    values.resize( numValues );
    prolongation.getRowCopy( globalRow, colIndices, values );

    for( localIndex i = 0; i < numValues; ++i )
    {
      globalIndex const colObj = colIndices[i] / numComp;
      integer const colComp = static_cast< integer >( colIndices[i] % numComp );
      views[colObj]( rowObj, colComp, rowComp ) = values[i];
    }
  }

  string_array fieldNames;
  fieldNames.insert( 0, names.begin(), names.end() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       fineManager,
                                                       domain.getNeighbors(),
                                                       false );

  writeFunc( names );
  for( string const & name : names )
  {
    fineManager.deregisterWrapper( name );
  }
}

}
}
}

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_
