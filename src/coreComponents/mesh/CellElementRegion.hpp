/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_CELLELEMENTREGION_HPP_
#define GEOSX_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

class CellElementRegionABC: public ElementRegionBaseABC
{
public:
  virtual localIndex getNElementsInCellElementSubRegion() const = 0;
  virtual std::list< std::reference_wrapper< const CellElementSubRegionABC > > getCellElementSubRegions() const = 0;
  virtual std::list< std::reference_wrapper< const CellElementSubRegionABC > > getElementSubRegionsSpec() const
  {
    return getCellElementSubRegions();
  }
  virtual localIndex getCellElementNbr() const = 0;
  virtual localIndex getSpecElementNbr() const
  {
    return getCellElementNbr();
  }
};

/**
 * @class CellElementRegion
 *
 * The CellElementRegion class contains the functionality to support the concept of a CellElementRegion in the element
 * hierarchy. CellElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class CellElementRegion : public ElementRegionBase, public CellElementRegionABC
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  CellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  CellElementRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~CellElementRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string catalogName()
  { return "CellElementRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return CellElementRegion::catalogName(); }

  ///@}

  /**
   * @name Generation of the cell element mesh region
   */
  ///@{

  /**
   * @brief Add a cellBlockRegion name to the list.
   * @param cellBlockName string containing the cell block region name.
   */
  void addCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.emplace_back( cellBlockName );
  }

  virtual void generateMesh( Group & cellBlocks ) override;

  /**
   * @brief Generate the aggregates.
   * @param faceManager a pointer to the FaceManager
   * @param nodeManager a pointer to the NodeManager
   */
  void generateAggregates( FaceManager const & faceManager, NodeManager const & nodeManager );

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the coarsening ratio
    static constexpr char const * coarseningRatioString() {return "coarseningRatio"; }

    /// @return String key for the cell block names
    static constexpr char const * sourceCellBlockNamesString() {return "cellBlocks"; }
  };

  /**
   * @brief Returns the number of elements in all the cell elememts sub-regions.
   * @return An integer
   *
   * @note virtual for testing reasons, we should have an ABC nevertheless.
   */
  virtual localIndex getNElementsInCellElementSubRegion() const override
  {
    return this->getNumberOfElements < CellElementSubRegion >();
  }

  /**
   * @brief Returns the CellElementSubRegion of the region manager.
   * @return An iterable to these elements.
   *
   * @note virtual for testing reasons, we should have an ABC nevertheless
   * @note Will become std::views with CXX20?
   */
  virtual std::list< std::reference_wrapper< const CellElementSubRegionABC > > getCellElementSubRegions() const override
  {
    std::list< std::reference_wrapper< const CellElementSubRegionABC > > result;

    auto append = [&result]( CellElementSubRegionABC const & er ) {
      result.push_back( er );
    };

    this->forElementSubRegions< CellElementSubRegionABC >( append );

    return result;
  }

private:

  // Cell block names
  string_array m_cellBlockNames;

  // Coarsening ratio
  real64 m_coarseningRatio;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTREGION_HPP_ */
