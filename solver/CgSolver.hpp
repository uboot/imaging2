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

#ifndef CGSOLVER_H
#define CGSOLVER_H

#include <solver/SolverInterface.hpp>

namespace imaging
{

  /** \ingroup solver
      \brief Wrapper class for the CG solver of provided by ITPACK.
      
      This class implements a conjugated gradients solver by interfacing <a href="http://rene.ma.utexas.edu/CNA/ITPACK">ITPACK</a>.
  */
  class CgSolver : public SolverInterface
  {
    size_t _n_max_iterations;

    void report_error(int error_flag) const;

  public:
    /** Constructs a CG solver. The solver iterates at most \em n_max_iterations times. */
    CgSolver(size_t n_max_iterations = 10000) : _n_max_iterations(n_max_iterations) {}

    void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const;
  };  
}

#endif
