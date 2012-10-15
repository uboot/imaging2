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

#ifndef SHAPE_MREP_XMLIO_H
#define SHAPE_MREP_XMLIO_H

#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>

#include <shape/mrep/Position2d.hpp>
#include <shape/mrep/MrepAtom.hpp>
#include <shape/mrep/MrepConnection.hpp>
#include <shape/mrep/MrepSkeleton2d.hpp>
#include <shape/mrep/MrepModel2d.hpp>
#include <shape/mrep/PolygonModel2d.hpp>

namespace imaging
{
  /** \ingroup mrep
      <tt>\#include <shape/mrep/xmlio.hpp></tt>
  */
  template<>
  class xml_handler< Position2d >
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, Position2d & object) const;
    void write_object(const Position2d & object, XmlWriter & out) const;

  };

  /** \ingroup mrep
      <tt>\#include <shape/mrep/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<MrepAtom>
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, MrepAtom & object) const;
    void write_object(const MrepAtom & object, XmlWriter & out) const;
  };

  /** \ingroup mrep
      <tt>\#include <shape/mrep/xmlio.hpp></tt>
  */
  template<>
  class xml_handler< MrepConnection<2> >
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, MrepConnection<2> & object) const;
    void write_object(const MrepConnection<2> & object, XmlWriter & out) const;
  };

  /** \cond */
  template<>
  class xml_handler<MrepSkeleton2d>
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, MrepSkeleton2d & object) const;
    void write_object(const MrepSkeleton2d & object, XmlWriter & out) const;
  };
  /** \endcond */

  /** \ingroup mrep
      <tt>\#include <shape/mrep/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<MrepModel2d>
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, MrepModel2d & object) const;
    void write_object(const MrepModel2d & object, XmlWriter & out) const;
  };

  /** \ingroup mrep
      <tt>\#include <shape/mrep/xmlio.hpp></tt>
  */
  template<>
  class xml_handler<PolygonModel2d>
  {
  public:
    static const std::string element_name;

    void read_object(XmlReader & in, PolygonModel2d & object) const;
    void write_object(const PolygonModel2d & object, XmlWriter & out) const;
  };

}


#endif
