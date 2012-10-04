#include <shape/xmlio.hpp>

namespace imaging
{
  const std::string xml_handler<Circle>::element_name = "circle";
    
  void xml_handler<Circle>::read_object(XmlReader & in, Circle & object) const
  {
    ublas::fixed_vector<float_t, 2> center;
    float_t radius;
    
    in >> XmlReader::attribute("center") >> center;
    in >> XmlReader::attribute("radius") >> radius;
    
    object.set_center(center);
    object.set_radius(radius);
  }
  
  void xml_handler<Circle>::write_object(const Circle & object, XmlWriter & out) const
  {
    out << XmlWriter::attribute("center") << object.center();
    out << XmlWriter::attribute("radius") << object.radius();
  } 
}

