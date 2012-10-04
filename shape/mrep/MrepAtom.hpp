// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MREPATOM_H
#define MREPATOM_H


namespace imaging
{
  /** \brief %MrepAtom */
  class MrepAtom
  {
    static const float_t RADIUS_COEFFICIENT = 5.0;
    float_t _radius;
    size_t _start_atom;
    size_t _connection;

  public:
    MrepAtom() : _radius(1.0), _start_atom(0), _connection(0) {}

    MrepAtom(float_t radius, size_t start_atom, size_t connection) :
        _radius(radius), _start_atom(start_atom), _connection(connection)
    {}

    void assign(float_t radius, size_t start_atom, size_t connection)
    {
      _radius = radius;
      _start_atom = start_atom;
      _connection = connection;

    }

    void set_radius(float_t radius) { _radius = radius; }
    void set_start_atom(size_t start_atom) { _start_atom = start_atom; }
    void set_connection(size_t connection) { _connection = connection; }
    
    float_t radius() const { return _radius; }
    size_t start_atom() const { return _start_atom; }
    size_t connection() const { return _connection; }

    void exponential(ublas::vector<float_t>::const_iterator & vector, MrepAtom & shape) const
    {
      shape._radius = _radius * exp(*vector); ++vector;
      // shape._radius = _radius + *vector; ++vector;
      shape._start_atom = _start_atom;
      shape._connection = _connection;
    }

    void logarithm(const MrepAtom & shape, ublas::vector<float_t>::iterator & vector) const
    {
      *vector = log( shape._radius /_radius ); ++vector;
      // *vector = shape._radius - _radius; ++vector;
    }

    size_t dimension() const { return 1; }
  };

  //   template <class T>
  //   std::ostream& operator<<(std::ostream & os, const MRepAtomS<T> & m_rep)
  //   {
  //     os << "MRepAtomS";
  //
  //     os << "\n";
  //
  //     os << "  radius: " << m_rep._radius << "\n";
  //
  //     os << std::endl;
  //
  //     return os;
  //   }

}


#endif



