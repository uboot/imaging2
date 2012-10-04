// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
