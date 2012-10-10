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

#ifndef POLYTOPE_SIMPLEPOLYGON_H
#define POLYTOPE_SIMPLEPOLYGON_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup polytope
      \brief Planar polygon defined by one polygonal boundary curve.
      
      This class implements a simple polygon given by the list of its boundary vertices. For more complex polygons consisting of multiple contours and holes use Polygon.
  */
  class SimplePolygon
  {
    std::auto_ptr< std::vector< ublas::fixed_vector<float_t, 2> > > _vertices;
    
    public:
    typedef std::vector< ublas::fixed_vector<float_t, 2> > vertex_list_t;
    
    SimplePolygon();
    
    /** Copy constructor. */
    SimplePolygon(const SimplePolygon & source);
    
    /** Copy assignement. */
    SimplePolygon & operator=(const SimplePolygon & source);
    
    /** Returns the <em>i</em>-th vertex of the polygon. */
    const ublas::fixed_vector<float_t, 2> & vertex(size_t i) const;
    
    /** Returns the total number of vertices of the polygon. */
    size_t n_vertices() const;
    
    /** Sets \em vertices as the vertices of the polygon. */
    void set_vertices(std::auto_ptr<vertex_list_t> vertices);
    
    /** Returns the vertices of the polygon. */
    const vertex_list_t & vertices() const;
  };
  
  /** \ingroup polytope 
      <tt>\#include <polytope/SimplePolygon.hpp></tt> 
      
      Computes the volume of the intersection of \em poly_a and \em poly_b. */
  float_t compute_intersection_volume(const SimplePolygon & poly_a, const SimplePolygon & poly_b);
}

#endif
