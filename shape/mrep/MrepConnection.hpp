// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MREPCONNECTION_H
#define MREPCONNECTION_H



namespace imaging
{
  /** \brief %MrepConnection */
  template<size_t N>
  class MrepConnection
    { };

  /** \brief %MrepConnection */
  template<>
  class MrepConnection<2>
  {
    float_t _rotation;
    float_t _radius;

  public:
    static const size_t DIMENSION = 2;

    MrepConnection<2>() : _rotation(0.0), _radius(1.0) {}

    MrepConnection<2>(float_t rotation, float_t radius) : _rotation(rotation), _radius(radius) {}

    void assign(float_t rotation, float_t radius) { _rotation = rotation, _radius = radius; }

    float_t radius() const { return _radius; }
    float_t rotation() const { return _rotation; }
    
    void set_radius(float_t radius) { _radius = radius; }
    void set_rotation(float_t rotation) { _rotation = rotation; }

    void exponential(ublas::vector<float_t>::const_iterator & vector, MrepConnection<2> & shape) const
    {
      shape._rotation = _rotation + *vector; ++vector;
      exponential_without_rotation(vector, shape);
    }

    void logarithm(const MrepConnection<2> & shape, ublas::vector<float_t>::iterator & vector) const
    {
      *vector = shape._rotation - _rotation; ++vector;
      logarithm_without_rotation(shape, vector);
    }

    void exponential_without_rotation(ublas::vector<float_t>::const_iterator & vector, MrepConnection<2> & shape) const
    {
      shape._radius = _radius * exp(*vector); ++vector;
      // shape._radius = _radius + *vector; ++vector;
    }

    void logarithm_without_rotation(const MrepConnection<2> & shape, ublas::vector<float_t>::iterator & vector) const
    {
      *vector = log( shape._radius /_radius ); ++vector;
      // *vector = shape._radius - _radius; ++vector;
    }

    size_t dimension() const { return 2; }
    size_t dimension_without_rotation() const { return 1; }
  };

  //   template <class S>
  //   std::ostream& operator<<(std::ostream & os, const MRepConnection<S> & connection)
  //   {
  //     os << "MRepConnection";
  //     os << "\n";
  //
  //     os << "  rotation: " <<  connection.get_rotation()<< "\n";
  //     os << "  distance: " << connection.get_distance() << "\n";
  //     os << std::endl;
  //
  //     return os;
  //   }
}



#endif
