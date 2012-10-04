#include <shape/mrep/Position2d.hpp>


namespace imaging
{
  Position2d::Position2d() :
  _center(0.0, 0.0), _rotation(0.0) {}

  Position2d::Position2d(const ublas::fixed_vector<float_t, 2> & center, float_t rotation) :
  _center(center), _rotation(rotation) {}

  void Position2d::assign(const ublas::fixed_vector<float_t, 2> & center, float_t rotation) { _center = center; _rotation = rotation; }

  void Position2d::exponential(ublas::vector<float_t>::const_iterator & vector, Position2d & shape) const
  {
    shape._center(0) = _center(0) + *vector; ++vector;
    shape._center(1) = _center(1) + *vector; ++vector;

    shape._rotation = _rotation + *vector; ++vector;
  }

  void Position2d::logarithm(const Position2d & shape, ublas::vector<float_t>::iterator & vector) const
  {
    *vector = shape._center(0) - _center(0); ++vector;
    *vector = shape._center(1) - _center(1); ++vector;

    *vector = shape._rotation - _rotation; ++vector;
  }
}




