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

#ifndef SHAPE_MREP_POSITION2D_H
#define SHAPE_MREP_POSITION2D_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \brief Stores the coordinates and the rotation of shapes in the plane. 
  
      This class stores coordinates in the plane and a rotation angle. Thus, it can be used to describe the position of shapes in the plane.
      Note that the members exponential() and logarithm() are similar but not exactly the same as in the class ShapeInterface.
  */
  class Position2d
  {
    ublas::fixed_vector<float_t, 2> _center;
    float_t _rotation;

  public:
    static const size_t DIMENSION = 2;

    Position2d();

    /** Constructs a Position2d object from \em center and \em rotation. */
    Position2d(const ublas::fixed_vector<float_t, 2> & center, float_t rotation);

    /** Assigns \em center and \em rotation to the object. */
    void assign(const ublas::fixed_vector<float_t, 2> & center, float_t rotation);

    /** Returns the center of the position. */
    const ublas::fixed_vector<float_t, 2> & center() const { return _center; }
    
    /** Returns the rotation of the position. */
    float_t rotation() const { return _rotation; }
    
    /** Sets the center of the position. */
    void set_center(const ublas::fixed_vector<float_t, 2> & center) { _center = center; }
    
    /** Sets the rotation of the position. */
    void set_rotation(float_t rotation) { _rotation = rotation; }

    /** Computes the exponential of the data at the current position of the iterator \em vector and stores it in \em shape. */
    void exponential(ublas::vector<float_t>::const_iterator & vector, Position2d & shape) const;

    /** Computes the logarithm of \em shape and stores it at the current position of the iterator \em vector. */
    void logarithm(const Position2d & shape, ublas::vector<float_t>::iterator & vector) const;
    
    /** Returns the data dimension of a position, i.e. \em 3. */
    size_t dimension() const { return 3; }
  };

}


#endif



