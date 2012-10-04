// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
