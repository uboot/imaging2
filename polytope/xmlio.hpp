// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef POLYTOPE_XMLIO_H
#define POLYTOPE_XMLIO_H

#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>
#include <polytope/Polygon.hpp>
#include <polytope/SimplePolygon.hpp>

namespace imaging
{
  /** \ingroup polytope
      <tt>\#include <polytope/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<Polygon>
  {
  public:
    /** \cond */
    static const std::string element_name;
    
    void read_object(XmlReader & in, Polygon & object) const;
    
    void write_object(const Polygon & object, XmlWriter & out) const;
    /** \endcond */
  }; 
  
  /** \ingroup polytope
      <tt>\#include <polytope/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<SimplePolygon>
  {
  public:
    /** \cond */
    static const std::string element_name;
    
    void read_object(XmlReader & in, SimplePolygon & object) const;
    
    void write_object(const SimplePolygon & object, XmlWriter & out) const;
    /** \endcond */
  }; 
}


#endif
