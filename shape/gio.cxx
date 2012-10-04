#include <shape/gio.hpp>

#include <spline/gio.hpp>

namespace imaging
{
  GraphicsInterface & operator<<(GraphicsInterface & out, const Circle & circle)
  {
    out.circle(circle.center(), circle.radius());
    
    return out;
  }

  
  GraphicsInterface & operator<<(GraphicsInterface & out, const BsplineShape & spline_shape)
  {
    out << spline_shape.curve();

    return out;
  }
}
