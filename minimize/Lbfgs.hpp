// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_LBFGS_H
#define MINIMIZE_LBFGS_H

#include <minimize/DifferentiableEnergyInterface.hpp>
#include <minimize/MinimizerInterface.hpp>

#include <external/liblbfgs/lbfgs.h>

namespace imaging
{

  /** \ingroup minimize
      \brief Minimizes differentiable energies using the Limited Memory BFGS method (a Quasi-Newton method).
      
      Minimizes a given energy by means of a Newton type method with estimated Hessian. This class and its documentation is based on <a href="http://www.chokkan.org/software/liblbfgs">liblbfgs</a>, a <em>C library of Limited memory BFGS (L-BFGS)</em> by Jorge Nocedal and Naoaki Okazaki.
  */
  class Lbfgs : public MinimizerInterface
  {
    DifferentiableEnergyInterface & _energy;
    lbfgs_parameter_t _lbfgs_parameters;
    size_t _n_steps;
    bool _terminated;
    
    static lbfgsfloatval_t evaluate(void *instance,
      const lbfgsfloatval_t *x,
      lbfgsfloatval_t *g,
      const int n,
      const lbfgsfloatval_t step);
      
    static int progress(void *instance,
      const lbfgsfloatval_t *x,
      const lbfgsfloatval_t *g,
      const lbfgsfloatval_t fx,
      const lbfgsfloatval_t xnorm,
      const lbfgsfloatval_t gnorm,
      const lbfgsfloatval_t step,
      int n,
      int k,
      int ls);

  public:
    /** Construct a Lbfgs object to minimize \em energy.  To actually start the minimization the user must call minimize().
    
        \param epsilon This parameter determines the accuracy with which the solution is to be found. A minimization terminates when ||g|| < \em epsilon * max(1, ||x||), where ||.|| denotes the Euclidean (L2) norm. The default value is 1e-5.
        \param line_search_f_tolerance A parameter to control the accuracy of the line search routine. The default value is 1e-4. This parameter should be greater than zero and smaller than 0.5. 
        \param line_search_g_tolerance The default value is 0.9. If the function gradient evaluations are inexpensive with respect to the cost of the iteration (which is sometimes the case when solving very large problems) it may be advantageous to set this parameter to a small value. A typical small value is 0.1. This parameter should be greater than the \em line_search_f_tolerance parameter (1e-4) and smaller than 1.0.
        \param n_correction_steps The number of corrections to approximate the inverse hessian matrix. The L-BFGS routine stores the computation results of previous \em n_correction_steps iterations to approximate the inverse hessian matrix of the current iteration. This parameter controls the size of the limited memories (corrections). The default value is 6. Values less than 3 are not recommended. Large values will result in excessive computing time.
        \param n_max_line_search_steps The maximum number of trials for the line search. This parameter controls the number of function and gradients evaluations per iteration for the line search routine. The default value is 20.  */
    Lbfgs(DifferentiableEnergyInterface & energy,
          float_t epsilon = 1.0e-5,
          float_t line_search_f_tolerance = 1.0e-4,
          float_t line_search_g_tolerance = 0.9,
          size_t n_correction_steps = 6,
          size_t n_max_line_search_steps = 20);

    bool minimize(size_t n_max_steps, size_t & n_steps); 
  };

}

#endif
