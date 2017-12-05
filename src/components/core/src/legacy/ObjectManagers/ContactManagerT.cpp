//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file ContactManagerT.cpp
 * @author Randolph Settgast
 * @date created on Sep 17, 2010
 */

#include "ContactManagerBaseT.h"
#include "ContactManagerT.h"

/**
 * @brief Constructor for the contacts / neighbor pairs
 * @author Scott Johnson
 * @date Jun 15, 2011
 * We are currently assuming that "contact" is the union of all "neighbor pairs"
 * and "inactive + active physical contacts"
 * Note that neighbor pairs that are also contacts would be the intersection of
 * these two sets
 * This data structure could be split into different structures later if it is
 * found that these are too cumbersome
 */
ContactManagerT::ContactManagerT():
  ContactManagerBaseT(ObjectDataStructureBaseT::ContactBaseManager),
  m_contactToIntersectionPolygonPointsMap(m_VariableOneToManyMaps["contactToIntersectionPolygonPointsMap"]),
  m_intersectionPolygonPoints(),
  ctr(0.0),
  dd(1e-4)
{
  this->AddKeylessDataField<R1Tensor>( "face1ParentSoln", true, true);
  this->AddKeylessDataField<R1Tensor>( "face2ParentSoln", true, true);
}

ContactManagerT::~ContactManagerT()
{}

localIndex ContactManagerT::NumberOfPolygons() const
{
  array<lArray1d>& tmp = m_contactToIntersectionPolygonPointsMap;
  localIndex num = 0;
  for( array<lArray1d>::size_type i = 0 ; i < tmp.size() ; ++i)
    num += tmp[i].size() > 0 ? 1 : 0;
  return num;
}

/**
 * @brief FOR DEFAULT VISUALIZATION - visit needs at least one polygon to
 * visualize a mesh or it gets angry
 * @author Scott Johnson
 * Calculates the center as the center of the domain and assigns the polygon a
 * dimension of 1e-4 the minimum dimension
 * @param[in] xmin Domain minima
 * @param[in] xmax Domain maxima
 */
void ContactManagerT::SetDefaultPolygonDimensions(const R1Tensor& xmin, const R1Tensor& xmax)
{
  R1Tensor dx = xmax;
  if(xmax.MaxVal() > 1e10 && xmin.MinVal() > 1e10)
  {
    this->ctr = 0.0;
    dx = 1e10;
  }
  else
  {
    this->ctr = xmax;
    this->ctr += xmin;
    this->ctr *= 0.5;
    dx -= xmin;
  }
  this->dd = dx.MinVal();
  this->dd *= 1e-4;
}
