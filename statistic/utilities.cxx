#include <statistic/utilities.hpp>

#include <core/utilities.hpp>

namespace imaging
{
  void mean(const ublas::matrix<float_t> & data, ublas::vector<float_t> & result)
  {
    result.resize(data.size2());
    result = ublas::scalar_vector<float_t>(result.size(), 0.0);
    
    for(ublas::matrix<float_t>::const_iterator1 iter1 = data.begin1(); iter1 != data.end1(); ++iter1)
      for(ublas::matrix<float_t>::const_iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        result(iter2.index2()) += *iter2;
        
    result = element_div(result, ublas::scalar_vector<float_t>(result.size(), float_t(data.size1())));
  }
  
  void var(const ublas::matrix<float_t> & data, ublas::vector<float_t> & result)
  {
    ublas::vector<float_t> data_mean;
    mean(data, data_mean);
    
    result.resize(data.size2());
    result = ublas::scalar_vector<float_t>(result.size(), 0.0);
    
    if(data.size1() == 1)
    {
      result(0) = 0.0;
      return;
    }
    
    for(ublas::matrix<float_t>::const_iterator1 iter1 = data.begin1(); iter1 != data.end1(); ++iter1)
      for(ublas::matrix<float_t>::const_iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        result(iter2.index2()) += square(*iter2 - data_mean(iter2.index2()));
        
    result = element_div(result, ublas::scalar_vector<float_t>(result.size(), float_t(data.size1() - 1)));
  }
  
  float_t mean(const ublas::vector<float_t> & data)
  {
    if(data.size() == 1)
      return data(0);
      
    float_t result = 0;
    
    for(ublas::vector<float_t>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
      result += *iter;
        
    result /= float_t(data.size());
    
    return result;
  }
  
  float_t var(const ublas::vector<float_t> & data)
  {
    if(data.size() == 1)
      return 0.0;
      
    float_t result = 0.0;
    float_t data_mean = mean(data);
    
    for(ublas::vector<float_t>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
      result += square(*iter - data_mean);
      
    result /= float_t(data.size());
      
    return result;
  }
}
