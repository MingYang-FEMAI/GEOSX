/*
 * LinearSolverWrapper.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#include "LinearSolverWrapper.hpp"

namespace geosx
{
namespace systemSolverInterface
{

LinearSolverWrapper::LinearSolverWrapper():
#if USE_MPI
  m_epetraComm(MPI_COMM_WORLD)
#else
  m_epetraComm()
#endif
{
  // TODO Auto-generated constructor stub

}

LinearSolverWrapper::~LinearSolverWrapper()
{
  // TODO Auto-generated destructor stub
}
}
} /* namespace geosx */
