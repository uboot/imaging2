// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_FUNCTIONALADAPTOR_H
#define MINIMIZE_FUNCTIONALADAPTOR_H

#include <core/FunctionalInterface.hpp>
#include <minimize/EnergyInterface.hpp>

namespace imaging
{
  template <class functional_t>
  class DifferentiableFunctionalAdaptor;
  
  /** \ingroup minimize 
      \brief Converts a FunctionalInterface object into an EnergyInterface object.
  */
  template <class functional_t>
  class FunctionalAdaptor : public EnergyInterface
  {
    friend class DifferentiableFunctionalAdaptor<functional_t>;
    
    functional_t & _functional;
    ublas::vector<float_t> _current_argument;
    float_t _current_energy;
    
  public:
    FunctionalAdaptor(functional_t & functional) :
      _functional(functional),
      _current_argument(functional.dimension()),
      _current_energy(0.0)
    {}
    
    ublas::vector<float_t> & current_argument()
    {
      return _current_argument;
    }
    
    void set_argument()
    {
      _current_energy = _functional(_current_argument);
    }
    
    float_t current_energy() const
    {
      return _current_energy;
    }
    
    std::size_t dimension() const
    {
      return _functional.dimension();
    }
  };
}

#endif
