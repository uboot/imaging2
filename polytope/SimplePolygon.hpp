// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
