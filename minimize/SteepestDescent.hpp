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

#ifndef MINIMIZE_STEEPESTDESCENT_H
#define MINIMIZE_STEEPESTDESCENT_H

#include <minimize/DifferentiableEnergyInterface.hpp>
#include <minimize/MinimizerInterface.hpp>

namespace imaging
{

  /** \ingroup minimize
      \brief Minimizes differentiable energies using steepest descent. 
      
      Minimizes a given energy \f$E\f$ which is initialized to a vector \f$x_0 \in \mathbf{R}^n\f$ by iteratively computing
      \f[
        x_{i+1} = -\tau \nabla E(x_i)\,,
      \f]
      for a step size \f$\tau > 0\f$.
  */
  class SteepestDescent : public MinimizerInterface
  {
    float_t _min_gradient_norm;
    float_t _step_size;
    bool _terminated;
    
    DifferentiableEnergyInterface & _energy;

  public:
    /** Construct a SteepestDescent object to minimize \em energy. If the L2-norm of the gradient of the functional falls below \em min_gradient_norm the algorithm stops. The parameter \em step_size refers to the step size \f$\tau\f$. To actually start the minimization the user must call minimize(). */
    SteepestDescent(DifferentiableEnergyInterface & energy, float_t min_gradient_norm, float_t step_size);
    
    bool minimize(size_t n_max_steps, size_t & n_actual_steps);
  };

}

#endif
