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

#ifndef BSPLINE_H
#define BSPLINE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup spline
      \brief Class template for B-Splines of arbitrary data type.
      
      This class template implements a B-Spline of arbitrary spline order. For template parameter \em DATA_t vector space operations (scalar multiplication, addition) and explicit construction from a float value must be implemented. Instantiating this class for 2-dimensional vectors creates a B-Spline curve.
      
      \sa PeriodicBspline
  */
  template <class DATA_t>
  class Bspline
  {
    size_t _spline_order;
    ublas::vector<DATA_t> _coefficients;
    ublas::vector<float_t> _knots;

    size_t determine_interval(float_t t) const
    {
      size_t i = 0;

      while( true )
      {
        if(_knots.size() == i)
          break;

        if(t < _knots(i))
          break;

        i++;
      }
      i = i - 1;

      if(i + _spline_order - 1 >= _knots.size() - 1 )
        i = _knots.size() - _spline_order - 1;

      return i;
    }

    bool is_in_spline_support(float_t t) const
    {
      return ( t >= first_knot() && t <= last_knot() );
    }

  public:
    Bspline() : _knots(), _coefficients() {}

    /** Construct a B-Spline of \em spline_order with \em n_coefficients. The number of knots is then the number of coefficients plus the spline order. */
    Bspline(size_t spline_order, size_t n_coefficients) : _spline_order(spline_order), _coefficients(n_coefficients + 1), _knots(n_coefficients + spline_order + 1) {}
    
    /** Copy constructor. */
    Bspline(const Bspline & source) : _spline_order(source._spline_order), _knots(source._knots), _coefficients(source._coefficients) {}

    /** Copy assignement. */
    const Bspline & operator=(const Bspline & rhs)
    { _spline_order = rhs._spline_order; _knots = rhs._knots; _coefficients = rhs._coefficients; return *this; }

    /** Returns the spline order. */
    size_t spline_order() const { return _spline_order; }
    
    /** Returns the <em>i</em>-th knot. */
    float_t knot(size_t i) const { return _knots(i); }
    
    /** Returns the first knot. */
    float_t first_knot() const { return _knots(_spline_order - 1); }
    
    /** Returns the last knot. */
    float_t last_knot() const { return _knots(n_knots() - _spline_order); }
    
    /** Sets the <em>i</em>-th knot to \em value. */
    void set_knot(size_t i, float_t value)
    {
#ifdef DEBUG
      if(i >= n_knots())
        throw Exception("Knot index out of bounds in Bspline::set_knot().");
#endif

      _knots(i) = value;

      if(i == n_knots() - 1)
        _knots(n_knots()) = value;
    }

    /** Returns the <em>i</em>-th coefficient. */
    const DATA_t & coefficient(size_t i) const { return _coefficients[i]; }
    
    /** Sets the <em>i</em>-th coefficient to \em value. */
    void set_coefficient(size_t i, const DATA_t & value)
    {
#ifdef DEBUG
      if(i >= n_coefficients())
        throw Exception("Coefficient index out of bounds in Bspline::set_coefficient().");
#endif

      _coefficients(i) = value;
    }

    /** Sets the order of the B-Spline to \em spline_order and its number of coefficients to \em n_coefficients. */
    void resize(size_t spline_order, size_t n_coefficients)
    {
      _spline_order = spline_order;
      _coefficients.resize(n_coefficients + 1);
      _knots.resize(n_coefficients + _spline_order + 1);
      _coefficients(_coefficients.size() - 1) = DATA_t(0.0);
    }

    /** Returns the number of knots of the B-Spline.  This is the number of coefficients plus the spline order. */
    size_t n_knots() const { return _knots.size() - 1; }
    
    /** Returns the number of coefficients of the B-Spline. */
    size_t n_coefficients() const { return _coefficients.size() - 1; }

    /** Evaluates the B-Spline at \em t and computes its \em first_derivative at \em t. If \em t is outside the spline support, \em first_derivative is set to 0 and 0 is returned. */
    DATA_t operator()(float_t t, DATA_t & first_derivative) const
      { DATA_t temp; return operator()(t, first_derivative, temp); }
      
    /** Evaluates the B-Spline at \em t. If \em t is outside the spline support, 0 is returned. */
    DATA_t operator()(float_t t) const
      { DATA_t temp_1, temp_2; return operator()(t, temp_1, temp_2); }

    /** Initializes the knots of the B-Spline to be equally distributed between \em x_0 and \em x_1. */ 
    void regular_knots(float_t x_0, float_t x_1)
    {
      if(n_coefficients() < _spline_order)
        throw Exception("Too little coefficients for Bspline::init_regular_knots().");


      float_t x = x_0;
      size_t i;

      for(i = 0; i < _spline_order; ++i)
        set_knot(i, x_0);

      for( ; i < n_knots() - _spline_order; ++i)
      {
        x += (x_1 - x_0) / float_t(n_knots() - 2 * _spline_order + 1);
        set_knot(i, x);
      }

      for(; i < n_knots(); ++i)
        set_knot(i, x_1);
    }

    /** Evaluates the B-Spline at \em t and computes its \em first_derivative and \em second_derivative at \em t. If \em t is outside the spline support, \em first_derivative and \em second_derivative are set to 0 and 0 is returned. */
    DATA_t operator()(float_t t, DATA_t & first_derivative, DATA_t & second_derivative) const
    {
      size_t k = _spline_order;

      if ( ! is_in_spline_support(t) )
      {
        first_derivative = DATA_t(0.0);
        second_derivative = DATA_t(0.0);
        return DATA_t(0.0);
      }

      size_t i = determine_interval(t);

      if(t < _knots(_spline_order - 1))
      {
        first_derivative = DATA_t(0.0);
        second_derivative = DATA_t(0.0);
        return DATA_t(0.0);
      }

      size_t i_offset = i - k + 1;
      ublas::matrix<float_t> basis_splines(k + 1, k);

      typename ublas::matrix<float_t>::iterator1 iter1;
      typename ublas::matrix<float_t>::iterator2 iter2;

      for(iter1 = basis_splines.begin1(); iter1 != basis_splines.end1(); ++iter1)
        for(iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
          *iter2 = 0.0;

      basis_splines(i - i_offset, 0) = 1;

      for(size_t j = 1; j < k; ++j)
        for(size_t ii = i - j; ii <= i; ++ii)
        {
          float_t n_1, n_2;

          if(_knots(ii + j) == _knots(ii))
            n_1 = 0;
          else
            n_1 = (t - _knots(ii)) /
                  (_knots(ii + j) - _knots(ii));

          if(_knots(ii + j + 1) == _knots(ii + 1))
            n_2 = 1;
          else
            n_2 = (_knots(ii + j + 1) - t) /
                  (_knots(ii + j + 1) - _knots(ii + 1));

          basis_splines(ii - i_offset, j) =
            basis_splines(ii - i_offset, j - 1) * n_1 +
            basis_splines(ii + 1 - i_offset, j - 1) * n_2;
        }


      DATA_t value = DATA_t(0.0);
      for (size_t ii = i - k + 1; ii <= i; ++ii)
        value += basis_splines(ii - i_offset, k - 1) * _coefficients(ii);


      ublas::vector<DATA_t> coefficient_differences(k); // corresponds to A^(1) in (13) in de Boor's paper
      size_t offset = i - k + 2; // offset between the 0 index of coefficient_differences and its real index.


      first_derivative = DATA_t(0.0);
      for (size_t ii = i - k + 2; ii <= i; ++ii)
      {
        if(_knots(ii + k - 1) != _knots(ii))
          coefficient_differences(ii - offset) = (_coefficients(ii) - _coefficients(ii - 1) ) /
                                                 (_knots(ii + k - 1) - _knots(ii) );
        else
          coefficient_differences(ii - offset) = DATA_t(0.0);
        first_derivative += basis_splines(ii - i_offset, k - 2)
                            * coefficient_differences(ii - offset);
      }
      first_derivative *= k - 1;

      second_derivative = DATA_t(0.0);
      for (size_t ii = i - k + 3; ii <= i; ++ii)
        if(_knots(ii + k - 2) != _knots(ii))
          second_derivative +=
            basis_splines(ii - i_offset, k - 3) /
            (_knots(ii + k - 2) - _knots(ii) )
            * (coefficient_differences(ii - offset) - coefficient_differences(ii - 1 - offset));
      second_derivative *= (k - 1) * (k - 2);


      return value;

    }
    
    /** Returns true if \em t is in the support of the <em>i</em>-th basis spline defined by the current knots of this spline. */
    bool is_in_basis_spline_support(size_t i, float_t t) const
    {
      return t >= knot(i) && t <= knot(i + spline_order());
    }
  };
}

#endif
