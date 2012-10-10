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
