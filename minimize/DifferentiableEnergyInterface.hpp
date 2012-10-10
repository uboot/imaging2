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

#ifndef MINIMIZE_DIFFERENTIABLEENERGYINTERFACE_H
#define MINIMIZE_DIFFERENTIABLEENERGYINTERFACE_H

#include <minimize/EnergyInterface.hpp>

namespace imaging
{
  /** \ingroup minimize 
      \brief Abstract class interface for differentiable energies depending on vector valued input data.
      
      In addition to the requirements of EnergyInterface, classes implementing DifferentiableEnergyInterface provide a current gradient computed from the current argument. These energy can be minimized using gradient based minimization techniques. 
  */
  class DifferentiableEnergyInterface : public EnergyInterface
  {
  public:
    /** Returns the current gradient. This function should \em not actually compute the current energy but return the cached result of the last call to set_argument()! */ 
    virtual const ublas::vector<float_t> & current_gradient() const = 0;
    
    /** Compute the energy value and the gradient corresponding to the current argument. This is where the main work of the energy evaluation should be done. In contrast to set_argument() this function also computes the gradient. I.e. in general it is more expensive than set_argument(). */
    virtual void set_argument_with_gradient() = 0;
  };
}

#endif
