// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_XMLIO_H
#define SHAPE_XMLIO_H

#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>

#include <shape/Circle.hpp>

namespace imaging
{ 
  /** \ingroup shape
      <tt>\#include <shape/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<Circle>
  {
  public:
    static const std::string element_name;
    
    void read_object(XmlReader & in, Circle & object) const;
    
    void write_object(const Circle & object, XmlWriter & out) const;
  }; 
}


#endif
