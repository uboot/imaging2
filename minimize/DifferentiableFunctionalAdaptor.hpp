// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MINIMIZE_DIFFERENTIABLEFUNCTIONALADAPTOR_H
#define MINIMIZE_DIFFERENTIABLEFUNCTIONALADAPTOR_H

#include <minimize/FunctionalAdaptor.hpp>
#include <minimize/DifferentiableEnergyInterface.hpp>

namespace imaging
{
  /** \ingroup minimize 
      \brief Converts a DifferentiableFunctionalInterface object into a DifferentiableEnergyInterface object.
  */
  template <class functional_t>
  class DifferentiableFunctionalAdaptor : public DifferentiableEnergyInterface
  {
    FunctionalAdaptor<functional_t> _functional_adaptor;
    ublas::vector<float_t> _current_gradient;
    
  public:
    DifferentiableFunctionalAdaptor(functional_t & functional) :
      _functional_adaptor(functional),
      _current_gradient(functional.dimension())
    {}
    
    ublas::vector<float_t> & current_argument()
    {
      return _functional_adaptor.current_argument();
    }
    
    void set_argument()
    {
      _functional_adaptor.set_argument();
    }
    
    void set_argument_with_gradient()
    {
      _functional_adaptor._current_energy = _functional_adaptor._functional(_functional_adaptor._current_argument, _current_gradient);
    }
    
    float_t current_energy() const
    {
      return _functional_adaptor.current_energy();
    }
    
    std::size_t dimension() const
    {
      return _functional_adaptor.dimension();
    }
    
    const ublas::vector<float_t> & current_gradient() const
    {
      return _current_gradient;
    }
  };
}

#endif
