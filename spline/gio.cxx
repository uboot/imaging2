#include <spline/gio.hpp>

namespace imaging
{
  GraphicsInterface & operator<<(GraphicsInterface & out, const Bspline< ublas::fixed_vector<float_t, 2> > & curve)
  {
    const size_t N_POINTS = 100;

    std::vector< ublas::fixed_vector<float_t, 2> > vertices(N_POINTS);

    float_t t = curve.first_knot();
    float_t offset = ( curve.last_knot() - curve.first_knot() ) / float_t(N_POINTS - 1);

    for(size_t i = 0; i < N_POINTS - 1; ++i)
    {
      vertices[i] = curve(t);
      t += offset;
    }
    
    vertices[N_POINTS - 1] = curve(curve.last_knot());

    out.polyline(vertices);

//     out.spline_curve(curve);

//     // draw control polygon
//     vertices.resize(curve.n_coefficients());
//     for(size_t i = 0; i < curve.n_coefficients(); ++i)
//       vertices[i] = curve.coefficient(i);
//     out.polygon(vertices);

    
    return out;
  }

  GraphicsInterface & operator<<(GraphicsInterface & out, const PeriodicBspline< ublas::fixed_vector<float_t, 2> > & curve)
  {
    const size_t N_POINTS = 100;

    std::vector< ublas::fixed_vector<float_t, 2> > vertices(N_POINTS);

    float_t t = curve.first_knot();
    float_t offset = ( curve.last_knot() - curve.first_knot() ) / float_t(N_POINTS);

    for(size_t i = 0; i < N_POINTS; ++i)
    {
      vertices[i] = curve(t);
      t += offset;
    }

    out.polygon(vertices);
    
//     // draw control polygon
//     vertices.resize(curve.n_coefficients());
//     for(size_t i = 0; i < curve.n_coefficients(); ++i)
//       vertices[i] = curve.coefficient(i);
//     out.polygon(vertices);

    return out;
  }

}

