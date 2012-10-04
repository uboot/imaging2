#include <polytope/xmlio.hpp>

#include <core/xmlio.hpp>

namespace imaging
{
  
  const std::string xml_handler<SimplePolygon>::element_name = "simple_polygon";
  const std::string xml_handler<Polygon>::element_name = "polygon";
  
  void xml_handler<SimplePolygon>::read_object(XmlReader & in, SimplePolygon & object) const
  {
    std::auto_ptr<SimplePolygon::vertex_list_t> vertices(new SimplePolygon::vertex_list_t);
    
    in >> XmlReader::element("vertex");
    in >> *vertices;
    
    object.set_vertices(vertices);
  }
    
  void xml_handler<SimplePolygon>::write_object(const SimplePolygon & object, XmlWriter & out) const
  {}
  
  
  void xml_handler<Polygon>::read_object(XmlReader & in, Polygon & object) const
  {
    std::auto_ptr< std::vector<SimplePolygon> > contours(new std::vector<SimplePolygon>);
    
    in >> XmlReader::element("contours");
    in >> *contours;
    in >> XmlReader::end_element;
    
//     in >> XmlReader::element("holes");
//     in >> object.holes();
//     in >> XmlReader::end_element;

    object.set_contours(contours);
  }
    
  void xml_handler<Polygon>::write_object(const Polygon & object, XmlWriter & out) const
  {}
}

