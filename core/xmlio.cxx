#include <core/xmlio.hpp>

#include <iostream>
#include <sstream>

#include <core/cio.hpp>

namespace imaging
{   
  const std::string xml_handler<ublas::matrix<float_t> >::element_name = "matrix";

  void xml_handler<ublas::matrix<float_t> >::read_object(XmlReader & in, ublas::matrix<float_t> & object) const
  {
    std::string data;
    in >> data;
    std::istringstream str_stream(data);
    str_stream >> object;
  }
  
  void xml_handler<ublas::matrix<float_t> >::write_object(const ublas::matrix<float_t> & object, XmlWriter & out) const
  {
    std::ostringstream str_stream;
    str_stream << object;
    out << str_stream.str();
  }

  const std::string xml_handler<ublas::vector<float_t> >::element_name = "vector";

  void xml_handler<ublas::vector<float_t> >::read_object(XmlReader & in, ublas::vector<float_t> & object) const
  {
    std::string data;
    in >> data;
    std::istringstream str_stream(data);
    str_stream >> object;
  }
  
  void xml_handler<ublas::vector<float_t> >::write_object(const ublas::vector<float_t> & object, XmlWriter & out) const
  {
    std::ostringstream str_stream;
    str_stream << object;
    out << str_stream.str();
  }
}

