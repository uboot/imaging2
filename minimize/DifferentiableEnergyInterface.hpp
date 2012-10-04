// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
