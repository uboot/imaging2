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

#ifndef POLYTOPE_POLYGON_H
#define POLYTOPE_POLYGON_H

#include <polytope/SimplePolygon.hpp>

namespace imaging
{

  /** \ingroup polytope
      \brief Planar polygon with multiple outlines and holes.
      
      This class implements a polygon defined by multiple contours and holes. Each of the contours and the holes is a SimplePolygon. For polygons consisting of one boundary contour only, use SimplePolygon instead.
  */
  class Polygon
  {
    std::auto_ptr< std::vector<SimplePolygon> > _contours;
    std::auto_ptr< std::vector<SimplePolygon> > _holes;
    
  public:
    Polygon();
    
    /** Copy constructor. */
    Polygon(const Polygon & source);
    
    /** Constructs a polygon object from a SimplePolygon. */
    explicit Polygon(const SimplePolygon & simple_polygon);
    
    /** Constructs a polygon from a vector of contours and a vector holes. */
    Polygon(std::auto_ptr< std::vector<SimplePolygon> > contours, std::auto_ptr< std::vector<SimplePolygon> > & holes);
    
    /** Copy assignement. */
    Polygon & operator=(const Polygon & source);
    
    /** Sets the contours of the polygon. */
    void set_contours(std::auto_ptr< std::vector<SimplePolygon> > contours);
    
    /** Sets the holes of the polygon. */
    void set_holes(std::auto_ptr< std::vector<SimplePolygon> > holes);
    
    /** Returns the contours of the polygon. */
    const std::vector<SimplePolygon> & contours() const { return *_contours; }
    
    /** Returns the holes of the polygon. */
    const std::vector<SimplePolygon> & holes() const { return *_holes; }
    
    /** Returns the number of contours of the polygon. */
    size_t n_contours() const { return _contours->size(); }
    
    /** Returns the number of holes of the polygon. */
    size_t n_holes() const { return _holes->size(); }
    
    /** Returns the <em>i</em>-th contour of the polygon. */
    const SimplePolygon & contour(size_t i) const { return (*_contours)[i]; }
    
    /** Returns the <em>i</em>-th hole of the polygon. */
    const SimplePolygon & hole(size_t i) const { return (*_holes)[i]; }
    
    /** Removes all contours and holes of the polygon. */
    void clear();
  };
  
  /** \ingroup polytope 
      <tt>\#include <polytope/Polygon.hpp></tt> 
      
      Computes the union of \em poly_1 and \em poly_2 and writes it to \em result. */
  void polygon_union(const Polygon & poly_1, const Polygon & poly_2, Polygon & result);

}

#endif
