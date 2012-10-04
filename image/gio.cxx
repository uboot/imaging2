#include <image/gio.hpp>

namespace imaging
{
  GraphicsInterface & operator<<(GraphicsInterface & out, const ColorImage2d & image)
  {
    out.image(image, ublas::fixed_vector<float_t, 2>(0.0, image.size()(0)), ublas::fixed_vector<float_t, 2>(0.0, image.size()(1)));
    
    return out;
  }

}
