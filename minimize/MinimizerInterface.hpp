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
