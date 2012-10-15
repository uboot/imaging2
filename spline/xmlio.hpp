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

#ifndef SPLINE_XMLIO_H
#define SPLINE_XMLIO_H

#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>
#include <spline/PeriodicBspline.hpp>
#include <core/xmlio.hpp>

namespace imaging
{
  /** \ingroup spline
      <tt>\#include <spline/xmlio.hpp></tt>
      
      Read and write periodic B-Spline curves from XML files.
      The format of the XML object is as follows:
  \code
  <periodic_spline_curve>
    <spline_order>3</spline_order>
    
    <!-- optional,
         if no knots are provided, equidistant knots starting at 0
         are used (as below) -->
    <knot>0.0</knot>
    <knot>1.0</knot>
    <knot>2.0</knot>
    <knot>3.0</knot>
    <knot>4.0</knot>
    
    <!-- use the right string represention of DATA_t
         within the coefficient-tags -->
    <coefficient>(0, 0)</coefficient> 
    <coefficient>(1, 0)</coefficient>
    <coefficient>(1, 1)</coefficient>
    <coefficient>(0, 1)</coefficient>
  </periodic_spline_curve>
  \endcode
  */
  template<>
  template <class DATA_t>
  class xml_handler< PeriodicBspline<DATA_t> >
  {
  public:
    /** \cond */
    static const std::string element_name;
    
    void read_object(XmlReader & in, PeriodicBspline<DATA_t> & object) const
    {
      std::vector<float_t> knots;
      std::vector<DATA_t> coefficients;
      size_t spline_order;
            
      in >> XmlReader::element("knot") >> knots;
      in >> XmlReader::element("coefficient") >> coefficients;
      in >> XmlReader::element("spline_order") >> spline_order >> XmlReader::end_element;
      
      object.resize(spline_order, coefficients.size());
      
      if(knots.size() == 0)
        for(size_t i = 0; i < object.n_knots(); ++i)
          object.set_knot(i, float_t(i));
      
      if(knots.size() != 0 && knots.size() != object.n_knots())
        throw Exception("Wrong number of knots in xml_handler<PeriodicBspline>::read_object().");
      
      for(size_t i = 0; i < coefficients.size(); ++i)
        object.set_coefficient(i, coefficients[i]);
      
      for(size_t i = 0; i < knots.size(); ++i)
        object.set_knot(i, knots[i]);
    }
    
    void write_object(const PeriodicBspline<DATA_t> & object, XmlWriter & out) const
    {
      out << XmlWriter::element("spline_order") << object.spline_order() << XmlWriter::end_element;
      
      out << XmlWriter::element("spline_order") << object.spline_order() << XmlWriter::end_element;
      
      for(size_t i = 0; i < object.n_coefficients(); ++i)
        out << XmlWriter::element("coefficient") << object.coefficient(i) << XmlWriter::end_element;
      
      for(size_t i = 0; i < object.n_knots(); ++i)
        out << XmlWriter::element("knot") << object.knot(i) << XmlWriter::end_element;
    }
    /** \endcond */
  }; 
    
  /** \cond */
  template<>
  template <class DATA_t>
  const std::string xml_handler< PeriodicBspline<DATA_t> >::element_name = "periodic_spline_curve";
  /** \endcond */
  
  
  /** \ingroup spline
      <tt>\#include <spline/xmlio.hpp></tt>
      
      Read and write B-Spline curves from XML files.
      The format of the XML object is as follows:
  \code
  <spline_curve>
    <spline_order>3</spline_order>
    
    <!-- optional,
         if no knots are provided, equidistant knots starting at 0
         are used (as below) -->
    <knot>0.0</knot>
    <knot>1.0</knot>
    <knot>2.0</knot>
    <knot>3.0</knot>
    <knot>4.0</knot>
    <knot>5.0</knot>
    <knot>6.0</knot>
    
    <!-- use the right string represention of DATA_t
         within the coefficient-tags -->
    <coefficient>(0, 0)</coefficient> 
    <coefficient>(1, 0)</coefficient>
    <coefficient>(1, 1)</coefficient>
    <coefficient>(0, 1)</coefficient>
  </spline_curve>
  \endcode
  */
  template<>
  template <class DATA_t>
  class xml_handler< Bspline<DATA_t> >
  {
  public:
    /** \cond */
    static const std::string element_name;
    
    void read_object(XmlReader & in, Bspline<DATA_t> & object) const
    {
      std::vector<float_t> knots;
      std::vector<DATA_t> coefficients;
      size_t spline_order;
            
      in >> XmlReader::element("knot") >> knots;
      in >> XmlReader::element("coefficient") >> coefficients;
      in >> XmlReader::element("spline_order") >> spline_order >> XmlReader::end_element;
      
      object.resize(spline_order, coefficients.size());
      
      if(knots.size() == 0)
        for(size_t i = 0; i < object.n_knots(); ++i)
          object.set_knot(i, float_t(i));
      
      if(knots.size() != 0 && knots.size() != object.n_knots())
        throw Exception("Wrong number of knots in xml_handler<Bspline>::read_object().");
      
      for(size_t i = 0; i < coefficients.size(); ++i)
        object.set_coefficient(i, coefficients[i]);
      
      for(size_t i = 0; i < knots.size(); ++i)
        object.set_knot(i, knots[i]);
    }
    
    void write_object(const Bspline<DATA_t> & object, XmlWriter & out) const
    {
      out << XmlWriter::element("spline_order") << object.spline_order() << XmlWriter::end_element;
      
      out << XmlWriter::element("spline_order") << object.spline_order() << XmlWriter::end_element;
      
      for(size_t i = 0; i < object.n_coefficients(); ++i)
        out << XmlWriter::element("coefficient") << object.coefficient(i) << XmlWriter::end_element;
      
      for(size_t i = 0; i < object.n_knots(); ++i)
        out << XmlWriter::element("knot") << object.knot(i) << XmlWriter::end_element;
    }
    /** \endcond */
  }; 
    
  /** \cond */
  template<>
  template <class DATA_t>
  const std::string xml_handler< Bspline<DATA_t> >::element_name = "spline_curve";
  /** \endcond */
}


#endif
