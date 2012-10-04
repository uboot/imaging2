#include <core/distribution_utilities.hpp>

#include <stdlib.h>
#include <boost/random.hpp>

namespace imaging
{
  namespace distribution_utilities_impl
  {
      boost::minstd_rand random_number_generator(static_cast<unsigned int>(time(0x0)));
      
      boost::uniform_real<float_t> unit_interval(0.0, 1.0);
      boost::uniform_real<float_t> symmetric_interval(-1.0, 1.0);
      boost::normal_distribution<float_t> std_normal(0.0, 1.0);
      
      boost::variate_generator<boost::minstd_rand&, boost::uniform_real<float_t> >
      uniform_distribution(random_number_generator, unit_interval);
      boost::variate_generator<boost::minstd_rand&, boost::uniform_real<float_t> >
      symmetric_uniform_distribution(random_number_generator, symmetric_interval);
      boost::variate_generator<boost::minstd_rand&, boost::normal_distribution<float_t> >
      std_normal_distribution(random_number_generator, std_normal);
  }
      
      
  float_t uniform_distribution()
  {
    return distribution_utilities_impl::uniform_distribution();
  }
  
  float_t symmetric_uniform_distribution()
  {
    return distribution_utilities_impl::symmetric_uniform_distribution();
  }
  
  float_t normal_distribution()
  {
    return distribution_utilities_impl::std_normal_distribution();
  }
}
