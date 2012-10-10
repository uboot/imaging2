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

#ifndef MINIMIZE_NLCG_H
#define MINIMIZE_NLCG_H

#include <minimize/DifferentiableEnergyInterface.hpp>
#include <minimize/MinimizerInterface.hpp>

namespace imaging
{

  /** \ingroup minimize
      \brief Minimizes differentiable energies using nonlinear conjugated gradients.
      
      Minimizes a given energy by means of the nonlinear conjugated gradients algorithm. This class is based on <a href="http://www.ece.northwestern.edu/%7Enocedal/CG%2B.html">CG+</a>, a <em>Software for Large-scale Unconstrained Optimization</em> by Guanghui Liu, Jorge Nocedal and Richard Waltz.
  */
  class NlCg : public MinimizerInterface
  {
    DifferentiableEnergyInterface & _energy;
    
    float_t _convergence_constant;
    int _method;
    ublas::vector<double> _d;
    ublas::vector<double> _g_old;
    ublas::vector<double> _w;
    int _iflag;
    int _finish;
    bool _terminated;

  public:
  
    enum methods {
      FLETCHER_REEVES /** Fletcher-Reeves */,
      POLAK_RIBIERE /** Polak-Ribiere */,
      POSITIVE_POLAK_RIBIERE  /** positive Polak-Ribiere */
    };
      
    /** Construct a NlCg object to minimize \em energy. The algorithm stops if the L2-norm of the solution falls below \em convergence_constant. To actually start the minimization the user must call minimize(). */
    NlCg(DifferentiableEnergyInterface & energy, float_t convergence_constant, size_t method = FLETCHER_REEVES);
    
    bool minimize(size_t n_max_steps, size_t & n_steps); 
  };

}

#endif
