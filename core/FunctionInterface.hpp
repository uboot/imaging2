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

#ifndef CORE_FUNCTIONINTERFACE_H
#define CORE_FUNCTIONINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class template for functions.
       
      Objects implementing this class template map objects of type \em A to objects of type \em B.
  */
  template <class A, class V>
  class FunctionInterface
  {
  public:
    virtual ~FunctionInterface() {}
  
    /** Evaluates the function at \em argument. */
    virtual void evaluate(const A & argument, V & value) = 0;
  };
}

#endif
