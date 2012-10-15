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

#ifndef SHAPE_CIRCLE_H
#define SHAPE_CIRCLE_H

#include <shape/ShapeInterface.hpp>
#include <shape/DiscretizableShapeInterface.hpp>

namespace imaging
{
  /** \ingroup shape
      \brief Class for circles in the plane.
      
      This class implements circular shapes in the plane.
  */
  class Circle : public ShapeInterface, public DiscretizableShapeInterface<2>
  {
    ublas::fixed_vector<float_t, 2> _center;
    float_t _radius;

    class Discretizer;
    
  public:
    const static size_t SHAPE_DIMENSION = 2;    

    Circle() : _center(0.0, 0.0), _radius(1.0) {}
    
    /** Construct a circle from \em center and \em radius. */
    Circle(const ublas::fixed_vector<float_t, 2> & center, float_t radius) : _center(center), _radius(radius) {}
    
    /** Copy constructor. */
    Circle(const Circle & source) : _center(source._center), _radius(source._radius) {}

    /** Copy assignement. */
    const Circle & operator=(const Circle & source);
    
    virtual std::auto_ptr< BoundaryDiscretizer<2> > boundary_discretizer(size_t n_points) const;

    /** Sets the \em center of the circle. */
    void set_center(const ublas::fixed_vector<float_t, 2> & center);
    
    /** Sets the \em radius of the circle. */
    void set_radius(float_t radius);
    
    /** Returns the center. */
    const ublas::fixed_vector<float_t, 2> & center() const { return _center; }
    
    /** Returns the radius. */
    float_t radius() const { return _radius; }

    void exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const;

    void logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const;
    
    size_t dimension() const { return 3; }
  };
  
  /** \cond */
  class Circle::Discretizer : public BoundaryDiscretizer<2>
  {
    float_t _step_size;
    const Circle & _circle;

  public:
    Discretizer(const Circle & circle, size_t n_points);

    void evaluate (size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & point, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const;
  };
  /** \endcond */
}


#endif
