#include <shape/Circle.hpp>

namespace imaging
{
  const Circle & Circle::operator=(const Circle & source)
  {
    _center = source._center; _radius = source._radius; return *this;
  }

  void Circle::set_center(const ublas::fixed_vector<float_t, 2> & center)
  {
    _center = center;
  }
  
  void Circle::set_radius(float_t radius)
  {
    _radius = radius;
  }
  
    
  void Circle::exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const
  {
    Circle * shape_ptr = dynamic_cast< Circle *>(& shape);
    
    if(! shape_ptr)
      throw Exception("Exception: Wrong type of argument 'shape' in Circle::exponential().");

    shape_ptr->_center(0) = _center(0) + vector(0);
    shape_ptr->_center(1) = _center(1) + vector(1);

    shape_ptr->_radius = _radius * exp(vector(2));
  }

  void Circle::logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const
  {
    const Circle * shape_ptr = dynamic_cast<const Circle *>(& shape);
    
    if(! shape_ptr)
      throw Exception("Exception: Wrong type of argument 'shape' in Circle::exponential().");

    vector(0) = shape_ptr->_center(0) - _center(0);
    vector(1) = shape_ptr->_center(1) - _center(1);

    vector(2) = log( shape_ptr->_radius / _radius );
  }
    
  
  Circle::Discretizer::Discretizer(const Circle & circle, size_t n_points) : BoundaryDiscretizer<2>(n_points), _circle(circle)
  {
    _step_size = 2 * PI / float_t(_n_points);
  }

  void Circle::Discretizer::evaluate (size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & point, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const
  {
    point(0) = _circle.center()(0) + _circle.radius()* cos(i * _step_size);
    point(1) = _circle.center()(1) + _circle.radius()* sin(i * _step_size);

    normal = point - _circle.center();
    normal *= _step_size;
    
    curvature = 1.0 / _circle.radius();
  }
  
  std::auto_ptr< BoundaryDiscretizer<2> > Circle::boundary_discretizer(size_t n_points) const
  {
    return std::auto_ptr< BoundaryDiscretizer<2> >(new Discretizer(*this, n_points));
  }
}


