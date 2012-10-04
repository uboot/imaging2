#include <core/vector_utilities.hpp>

#include <core/utilities.hpp>

namespace imaging
{
  const ublas::fixed_vector<float_t, 2> polar2cartesian(float_t radius, float_t angle)
  {
    ublas::fixed_vector<float_t, 2> result;
    
    result(0) = radius * cos(angle);
    result(1) = radius * sin(angle);

    return result;
  }
  
  float_t radius(const ublas::fixed_vector<float_t, 2> & v)
  {
    return norm_2(v);
  }

  float_t angle(const ublas::fixed_vector<float_t, 2> & v)
  {
    float_t local_angle;
    float_t abs_value = radius(v);

    if(v(0) != 0.0)
      if(v(0) > 0.0)
        local_angle = asin(v(1) / abs_value);
      else
        local_angle = PI - asin(v(1) / abs_value);
    else
      if(v(1) != 0)
        if(v(1) > 0.0)
          local_angle = PI / 2.0;
        else
          local_angle = - PI / 2.0;
      else
        local_angle = 0.0;

    if(local_angle < 0) local_angle += 2 * PI;

    return local_angle;
  }
}

