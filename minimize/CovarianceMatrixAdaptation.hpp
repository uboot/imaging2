// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_COVARIANCEMATRIXADAPTATION_H
#define MINIMIZE_COVARIANCEMATRIXADAPTATION_H

#include <minimize/EnergyInterface.hpp>
#include <minimize/MinimizerInterface.hpp>

namespace imaging
{
  /** \ingroup minimize
      \brief Minimizes energies using covariance matrix adaptation (CMA).
      
      This class attempts (module programming errors) to implement the covariance matrix adaptation as in <em><a href="http://www.bionik.tu-berlin.de/user/niko/cmatutorial.pdf">Nikolaus Hansen, "The CMA Evolution Strategy: A Tutorial"</a></em>.
  */
  class CovarianceMatrixAdaptation : public MinimizerInterface
  {
    std::size_t _lambda;
    float_t _sigma;
    std::size_t _mu;
    ublas::vector<float_t> _w;
    float_t _c_sigma;
    float_t _d_sigma;
    float_t _mu_eff;
    float_t _mu_cov;
    float_t _c_cov;
    float_t _c_c;
    std::size_t _n;
    ublas::matrix<float_t> _C;
    ublas::vector<float_t> _p_sigma;     
    ublas::vector<float_t> _p_c; 
    bool _terminated;
      
    EnergyInterface & _energy;
    
    float_t _min_update;
    
    void init(float_t sigma, float_t min_update, std::size_t lambda);
  
  public:
    /** Construct a CMA object to minimize \em energy. The minimum should not be more 3 \em sigma away from the current argument of \em energy. The parameter \em lambda denotes the population size in the CMA algorithm. If the threshold \em minimal_update is not met during a step of the minimization the algorithm stops. To actually start the minimization the user must call minimize(). */
    CovarianceMatrixAdaptation(EnergyInterface & energy, float_t sigma, float_t min_update, std::size_t lambda);  
      
    /** Construct a CMA object to minimize \em energy. The minimum should not be more 3 \em sigma away from the current argument of \em energy. If the threshold \em minimal_update is not met during a step of the minimization the algorithm stops. To actually start the minimization the user must call minimize(). */
    CovarianceMatrixAdaptation(EnergyInterface & energy, float_t sigma, float_t min_update); 

    bool minimize(size_t n_max_steps, size_t & n_actual_steps);
  };

}

#endif
