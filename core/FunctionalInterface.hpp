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

#ifndef CORE_FUNCTIONALINTERFACE_H
#define CORE_FUNCTIONALINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class interface for functionals depending on vector valued input data.
  */
  class FunctionalInterface
  {
  public:
    virtual ~FunctionalInterface() {}
  
    /** Evaluates the functional at \em x. */
    virtual float_t operator()(const ublas::vector<float_t> & x) = 0;
    
    /** Returns the dimension of the input argument for this functional. In other words, this member returns the dimension of the space the functional is defined on. */
    virtual std::size_t dimension() const = 0;
  };
}

#endif
