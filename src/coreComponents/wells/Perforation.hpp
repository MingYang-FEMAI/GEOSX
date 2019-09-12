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

/*
 * @file Perforation.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_PERFORATION_HPP
#define GEOSX_CORECOMPONENTS_WELLS_PERFORATION_HPP

#include "dataRepository/Group.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforation = "Perforation";
}
}


/**
 * @class Perforation
 *
 * This class describes a perforation with its location, transmissibility and corresponding well element
 */  
class Perforation : public dataRepository::Group
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  explicit Perforation( string const & name, dataRepository::Group * const parent );
  
  /**
   * @brief default destructor
   */
  ~Perforation() override;

  /// deleted default constructor
  Perforation() = delete;

  /// deleted copy constructor
  Perforation( Perforation const &) = delete;

  /// deleted move constructor
  Perforation( Perforation && ) = delete;

  /// deleted assignment operator
  Perforation & operator=( Perforation const & ) = delete;

  /// deleted move operator
  Perforation & operator=( Perforation && ) = delete;

  /**
   * @brief Getter for the linear distance between the well head and the perforation
   * @return the distance between the well head and the perforation
   */
  real64 const & GetDistanceFromWellHead() const { return m_distanceFromHead; }

  /**
   * @brief Getter for the transmissibility at the perforation
   * @return the transmissibility
   */
  real64 GetTransmissibility() const { return m_transmissibility; }


  struct viewKeyStruct
  {
    static constexpr auto distanceFromHeadString = "distanceFromHead";
    static constexpr auto transmissibilityString = "transmissibility";

    dataRepository::ViewKey distanceFromHead = { distanceFromHeadString };
    dataRepository::ViewKey transmissibility = { transmissibilityString };

  } viewKeysPerforation;

protected:

  void PostProcessInput() override;

private:
  
  // linear distance from well head
  real64 m_distanceFromHead;

  // transmissibility (well index) at this perforation
  real64 m_transmissibility;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
