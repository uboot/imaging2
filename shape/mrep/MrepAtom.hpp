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



