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

#ifndef SOLVER_LUSOLVER_H
#define SOLVER_LUSOLVER_H

#include <solver/SolverInterface.hpp>
                    
namespace imaging
{                                      
  /** \ingroup solver
      \brief Wrapper class for the LU solver provided by SPOOLES.
      
      This class implements a conjugated gradients solver by interfacing <a href="http://www.netlib.org/linalg/spooles/spooles.2.2.html">SPOOLES</a>.
  */
  class LuSolver : public SolverInterface
  {

  public:
    LuSolver() {}
    
    void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const;
  };

}

#endif
