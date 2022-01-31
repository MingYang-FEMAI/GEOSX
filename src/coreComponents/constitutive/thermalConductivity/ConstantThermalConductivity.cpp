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
 * @file ConstantThermalConductivity.cpp
 */

#include "ConstantThermalConductivity.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ConstantThermalConductivity::ConstantThermalConductivity( string const & name, Group * const parent ):
  ThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::thermalConductivityComponentsString(), &m_thermalConductivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz components of a diagonal thermal conductivity tensor [W/(m.K)]" );
}

std::unique_ptr< ConstitutiveBase >
ConstantThermalConductivity::deliverClone( string const & name,
                                           Group * const parent ) const
{
  return ThermalConductivityBase::deliverClone( name, parent );
}

void ConstantThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      m_effectiveConductivity[ei][q][0] = m_thermalConductivityComponents[0];
      m_effectiveConductivity[ei][q][1] = m_thermalConductivityComponents[1];
      m_effectiveConductivity[ei][q][2] = m_thermalConductivityComponents[2];
    }
  }
}

void ConstantThermalConductivity::postProcessInput()
{
  GEOSX_THROW_IF( m_thermalConductivityComponents[0] <= 0 ||
                  m_thermalConductivityComponents[1] <= 0 ||
                  m_thermalConductivityComponents[2] <= 0,
                  GEOSX_FMT( "{}: the components of the thermal conductivity tensor must be strictly positive",
                             getFullName() ),
                  InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geosx