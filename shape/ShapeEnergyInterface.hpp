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

#ifndef SHAPE_SHAPEENERGYINTERFACE_H
#define SHAPE_SHAPEENERGYINTERFACE_H

#include <minimize/EnergyInterface.hpp>


namespace imaging
{
  /** \ingroup shape 
      \brief Abstract class interface for energies depending on vector valued input data which define an input shape.
      
      In addition to the requirements of EnergyInterface, classes implementing ShapeEnergyInterface provide a current shape computed from the current argument. Minimizers can access this shape to obtain geometric information (via BoundaryDiscretizer) about the current argument. 
  */
  template <class shape_t>
  class ShapeEnergyInterface : public EnergyInterface
  {
  public:
  
    /** Returns the current shape. This function should \em not actually compute the current shape but return the cached result of the last call to set_argument()! */ 
    virtual const shape_t & current_shape() const = 0;
  };
}

#endif
