// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
