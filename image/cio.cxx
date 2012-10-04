#include <image/cio.hpp>

#include <core/cio.hpp>
#include <iostream>

namespace imaging
{ 
  std::ostream & operator<<(std::ostream & out, const GrayValue value)
  {
    unsigned char temp(value);
    return out << (unsigned int)(temp);
  }
  
  std::ostream & operator<<(std::ostream & out, const Color & value)
  {
    return out << (ublas::fixed_vector<GrayValue, 3> &)(value);
  }

  std::istream & operator>>(std::istream & in, GrayValue & value)
  {
    unsigned int temp;
    in >> temp;
    value = (unsigned char &)(temp);
    return in;
  }
  
  std::istream & operator>>(std::istream & in, Color & value)
  {
    return in >> (ublas::fixed_vector<GrayValue, 3> &)(value);
  }
}
