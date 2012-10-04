#include <shape/BsplineShape.hpp>

namespace imaging
{
  void BsplineShape::exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const
  {
    BsplineShape * shape_ptr = dynamic_cast<BsplineShape *>(& shape);
    
    if(! shape_ptr)
      throw Exception("Exception: Wrong type of argument 'shape' in BsplineShape::exponential().");

    shape_ptr->_curve = _curve;
    for(size_t i = 0; i < _curve.n_coefficients(); ++i)
    {
      ublas::fixed_vector<float_t, 2> offset;
      offset(0) = vector(2*i);
      offset(1) = vector(2*i + 1);
      shape_ptr->_curve.set_coefficient(i, _curve.coefficient(i) + offset);
    }
  }

  void BsplineShape::logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const
  { 
    const BsplineShape * shape_ptr = dynamic_cast<const BsplineShape *>(& shape);
    
    if(! shape_ptr)
      throw Exception("Exception: Wrong type of argument 'shape' in BsplineShape::exponential().");

    for(size_t i = 0; i < _curve.n_coefficients(); ++i)
    {
      vector(2*i) = shape_ptr->_curve.coefficient(i)(0) - _curve.coefficient(i)(0);
      vector(2*i + 1) = shape_ptr->_curve.coefficient(i)(1) - _curve.coefficient(i)(1);
    }
  }

  BsplineShape::Discretizer::Discretizer(const BsplineShape  & spline_curve, size_t n_points) : BoundaryDiscretizer<2>(n_points), _spline_curve(spline_curve)
  {
    _step_size = (spline_curve._curve.last_knot() - spline_curve._curve.first_knot()) / float_t(_n_points);
  }

  void BsplineShape::Discretizer::evaluate (size_t i, ublas::fixed_vector<float_t, 2> & point, ublas::fixed_vector<float_t, 2> & normal, float_t & curvature) const
  { 
    ublas::fixed_vector<float_t, 2> tangent, second_derivative;
    point = _spline_curve._curve(i * _step_size, tangent, second_derivative);

    tangent *= _step_size;

    normal(0) = - tangent(1);
    normal(1) = tangent(0);
    
    curvature = boundary_discretizer_impl::compute_curve_curvature(tangent, second_derivative);
  }
  
  std::auto_ptr< BoundaryDiscretizer<2> > BsplineShape::boundary_discretizer(size_t n_points) const
  {
    return std::auto_ptr< BoundaryDiscretizer<2> >(new Discretizer(*this, n_points));
  }
}

