// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef PERIODICBSPLINE_H
#define PERIODICBSPLINE_H


#include <spline/Bspline.hpp>

namespace imaging
{
  /** \ingroup spline
      \brief Class template for periodic B-Splines of arbitrary data type.
      
      This class template implements a periodic B-Spline of spline arbitrary order. For template parameter \em DATA_t vector space operations (scalar multiplication, addition) and explicit construction from a float value must be implemented. Instantiating this class for 2-dimensional vectors creates a closed B-Spline curve.
      
      The last knot of a periodic B-spline is identified with the first one, i.e. the number of knots equals the number of coefficients plus 1. E.g. the knot vector [1, 2, 3, 4] corresponds to an equilateral "knot triangle".
      
      \sa Bspline
  */
  template <class DATA_t>
  class PeriodicBspline : Bspline<DATA_t>
  {
    float_t transform_parameter(float_t t) const
    {
      // scale to [0, 1]
      float_t transformed_t = ( t - first_knot() ) / ( last_knot() - first_knot() );
      
      // move to [0, 1]
      transformed_t = transformed_t - floor(transformed_t);
      
      //rescale to orginal domain
      return transformed_t * ( last_knot() - first_knot() ) + knot(0);
    }
    
  public:
    PeriodicBspline() : Bspline<DATA_t>() {}

    /** Copy constructor. */
    PeriodicBspline(const PeriodicBspline & source) : Bspline<DATA_t>(source) {}

    /** Construct a periodic B-Spline of \em spline_order with \em n_coefficients. The number of knots is then the number of coefficients plus 1. */
    PeriodicBspline(size_t spline_order, size_t n_coefficients)
    { resize(spline_order, n_coefficients); }

    /** Copy assignement. */
    const PeriodicBspline & operator=(const PeriodicBspline & rhs)
    { Bspline<DATA_t>::operator=(rhs); return *this; }

    size_t spline_order() const { return Bspline<DATA_t>::spline_order(); }
    float_t knot(size_t i) const { return Bspline<DATA_t>::knot(i + Bspline<DATA_t>::spline_order() - 1); }

    float_t first_knot() const { return Bspline<DATA_t>::first_knot(); }
    float_t last_knot() const { return Bspline<DATA_t>::last_knot(); }

    void set_knot(size_t i, float_t value)
    {
      Bspline<DATA_t>::set_knot(i + Bspline<DATA_t>::spline_order() - 1, value);

      for(size_t j = n_knots(); j < n_knots() + Bspline<DATA_t>::spline_order() - 1; ++j)
        Bspline<DATA_t>::set_knot(j + Bspline<DATA_t>::spline_order() - 1, Bspline<DATA_t>::knot(j - 1 + Bspline<DATA_t>::spline_order() - 1) +
                                        ( knot(j - n_knots() + 1) - knot(j - n_knots()) ) );

      for(size_t j = 1; j < Bspline<DATA_t>::spline_order(); ++j)
        Bspline<DATA_t>::set_knot(Bspline<DATA_t>::spline_order() - 1 - j, Bspline<DATA_t>::knot(Bspline<DATA_t>::spline_order() - 1 - j + 1) -
                                        ( knot(n_knots() - j) - knot(n_knots() - j - 1) ) );
    }

    void set_coefficient(size_t i, const DATA_t & value)
    {
      if( i + Bspline<DATA_t>::spline_order() > n_coefficients() )
        Bspline<DATA_t>::set_coefficient(i + Bspline<DATA_t>::spline_order() - 1 - n_coefficients(), value);

      Bspline<DATA_t>::set_coefficient(i + Bspline<DATA_t>::spline_order() - 1, value);
    }

    const DATA_t & coefficient(size_t i) const { return Bspline<DATA_t>::coefficient(i + Bspline<DATA_t>::spline_order() - 1); }

    void resize(size_t spline_order, size_t n_coefficients)
    {
      if(n_coefficients + 1 < spline_order)
        throw Exception("Exception: Too little coefficients in PeriodicBspline::resize().");

      Bspline<DATA_t>::resize(spline_order, n_coefficients + spline_order - 1);
    }

    size_t n_coefficients() const
    { return Bspline<DATA_t>::n_coefficients() - Bspline<DATA_t>::spline_order() + 1; }

    /** Returns the number of knots of the B-Spline.  This is the number of coefficients plus 1. */
    size_t n_knots() const
      { return n_coefficients() + 1; }

    DATA_t operator()(float_t t, DATA_t & first_derivative) const
      { DATA_t temp; return operator()(t, first_derivative, temp); }

    DATA_t operator()(float_t t) const
      { DATA_t temp_1, temp_2; return operator()(t, temp_1, temp_2); }

    DATA_t operator()(float_t t, DATA_t & first_derivative, DATA_t & second_derivative) const
    {
      float_t transformed_t = transform_parameter(t);

      DATA_t value = Bspline<DATA_t>::operator()(transformed_t, first_derivative, second_derivative);

      return value;
    }

    void regular_knots(float_t x_0, float_t x_1)
    {
      float_t offset = (x_1 - x_0) / float_t(n_coefficients());
      float_t x = x_0;
      set_knot(0, x);
      for(size_t i = 1; i < n_knots(); ++i)
      {
        x += offset;
        set_knot(i, x);
      }
    }
    
    /** Returns true if \em t is in the support of the <em>i</em>-th basis spline defined by the current knots of this spline. */
    bool is_in_basis_spline_support(size_t i, float_t t) const
    {
      float_t transformed_t = transform_parameter(t);
      size_t looped_right_knot_index = ( i + spline_order() ) % n_knots();
      return t >= max(knot(i), looped_right_knot_index) && t <= knot(i + spline_order());
    }
    
  };
}

#endif
