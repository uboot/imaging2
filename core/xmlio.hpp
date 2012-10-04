// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_XMLIO_H
#define CORE_XMLIO_H

#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>

namespace imaging
{ 
  /** \ingroup core
      <tt>\#include <core/xmlio.hpp></tt>
  */
  template<>
  class xml_handler< ublas::matrix<float_t> >
  {
  public:
    static const std::string element_name;
    
    void read_object(XmlReader & in, ublas::matrix<float_t> & object) const;
    
    void write_object(const ublas::matrix<float_t> & object, XmlWriter & out) const;
  }; 
  
  /** \ingroup core
      <tt>\#include <core/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<ublas::vector<float_t> >
  {
  public:
    static const std::string element_name;
    
    void read_object(XmlReader & in, ublas::vector<float_t> & object) const;
    
    void write_object(const ublas::vector<float_t> & object, XmlWriter & out) const;
  }; 
}


#endif
