// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_MINIMIZERINTERFACE_H
#define MINIMIZE_MINIMIZERINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup minimize 
      \brief Abstract base class of iterative minimization methods. 
      
      This interface is a common base class of iterative minimization methods. The parameters of a specific minimization method should be passed in the constructor together with a reference to the energy to be minimized. The actual iterative minimization is started by calling minimize(). It should be possible to call minimize() several times in case the convergence criterion was not met during the first call. However, due to the complexity of some minimization routines it is in general not possible to guarantee that iterative calls of minimize() behave exactly the same as one call with a higher value \em n_max_steps.
  */
  class MinimizerInterface
  {
  public:
    virtual ~MinimizerInterface() {};
    
      /** Start the minimization process. At most \em n_max_steps will be performed. If the convergence criterion (as set in the constructor) is met before the number of maximal steps is reached the function returns true. Otherwise it returns false. The actual number of performed stops is stored in \em n_steps. Upon return, the energy (as set in the constructor) will have the solution as its current argument. */
    virtual bool minimize(size_t n_max_steps, size_t & n_actual_steps) = 0;
  };
}


#endif
