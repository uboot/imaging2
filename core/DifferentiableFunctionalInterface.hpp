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

#ifndef CORE_DIFFERENTIABLEFUNCTIONALINTERFACE_H
#define CORE_DIFFERENTIABLEFUNCTIONALINTERFACE_H

#include <core/FunctionalInterface.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class interface of differentiable functionals depending on vector valued input data. */
  class DifferentiableFunctionalInterface : public FunctionalInterface
  {
  public:
    using FunctionalInterface::operator();
    
    virtual float_t operator()(const ublas::vector<float_t> & x, ublas::vector<float_t> & gradient) = 0;
  };
}

#endif
