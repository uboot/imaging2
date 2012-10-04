#ifndef CORE_GIO_H
#define CORE_GIO_H

#include <graphics/GraphicsInterface.hpp>

namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/gio.hpp></tt>
  */
  GraphicsInterface & operator<<(GraphicsInterface & out, const ublas::fixed_vector<float_t, 2> & vertex)
  {
    out.vertex(vertex);
    
    return out;
  }

}


#endif
