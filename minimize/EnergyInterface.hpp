// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_ENERGYINTERFACE_H
#define MINIMIZE_ENERGYINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup minimize 
      \brief Abstract class interface for energies depending on vector valued input data.
      
      Classes implementing EnergyInterface always have a current argument which can be directly accessed from the outside. It is only when set_argument() is called that the energy value for the current argument is computed. In other words, instead of passing an argument to a member of the energy class the clients edit the argument inside the energy class and manually trigger its evaluation. This has the advantage that the argument has only to be stored once and is still always accesible from the outside. Minor modifications of the argument can be made without resetting it completely. Accessing the argument is expected to be cheap, only calls to set_argument() might be expensive.
  */
  class EnergyInterface
  {
  public:
    virtual ~EnergyInterface() {}
  
    /** Access the current argument of the energy. In general the user should not resize the current argument! */
    virtual ublas::vector<float_t> & current_argument() = 0;
    
    /** Compute the energy value corresponding to the current argument. This is where the main work of the energy evaluation should be done. */
    virtual void set_argument() = 0;
    
    /** Returns the current energy. This function should \em not actually compute the current energy but return the cached result of the last call to set_argument()! */
    virtual float_t current_energy() const = 0;
    
    /** Returns the dimension the class expects as input data. The function current_argument() will return a vector of this dimension and the user should not change the size of this vector! */
    virtual std::size_t dimension() const = 0;
  };
}

#endif
