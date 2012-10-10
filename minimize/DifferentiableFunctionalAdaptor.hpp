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
