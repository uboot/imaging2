#include <polytope/gio.hpp>

namespace imaging
{
  GraphicsInterface & operator<<(GraphicsInterface & out, const Polygon & polygon)
  {
    for(size_t i = 0; i < polygon.n_holes(); ++i)
      out << polygon.hole(i);
      
    for(size_t i = 0; i < polygon.n_contours(); ++i)
      out << polygon.contour(i);
      
    return out;
  } 
  
  GraphicsInterface & operator<<(GraphicsInterface & out, const SimplePolygon & polygon)
  {
    out.polygon(polygon.vertices());
    
    return out;
  }
}

