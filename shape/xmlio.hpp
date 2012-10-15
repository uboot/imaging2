/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

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
