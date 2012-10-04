// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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



