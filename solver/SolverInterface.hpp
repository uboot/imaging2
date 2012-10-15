/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOLVER_SOLVERINTERFACE_H
#define SOLVER_SOLVERINTERFACE_H

#include <core/imaging2.hpp>
                    
namespace imaging
{                                      
  /** \ingroup solver
      \brief Abstract class interface for solver for sparse systems of linear equations.
  */
  class SolverInterface
  {

  public:
    virtual ~SolverInterface() {}
  
    /** Solves the system of equations defined by the matrix \em eqs and the vector \em rhs and writes the solution to \em result. The vector \em result is automatically resized to dimension of the system. If the sizes of the input data are such that the system cannot be solved an Exception is thrown. */
    virtual void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const = 0;
  };

}

#endif
