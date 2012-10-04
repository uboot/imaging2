#include <core/utilities.hpp>


namespace imaging
{
  // TODO: review function!!!
  float_t clockwise_difference(float_t angle_1, float_t angle_2)
  {
    if(angle_1 >= angle_2 - 0.0001)
      return angle_1 - angle_2;
    else
      return angle_1 + (2 * PI - angle_2);
  }
  
  float_t counter_clockwise_difference(float_t angle_1, float_t angle_2)
  {
    return clockwise_difference(angle_2, angle_1);
  }
  
  float_t delta(size_t base, size_t t)
  {
    if(base == t)
      return 1.0;
    else 
      return 0.0;
  }

}
