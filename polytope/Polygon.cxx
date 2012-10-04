#include <polytope/Polygon.hpp>

extern "C"
{
#include <external/gpc/gpc.h>
}

namespace imaging
{
  Polygon::Polygon()
  {
    _contours = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>);
    _holes = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>);
  }
  
  Polygon::Polygon(const Polygon & source)
  {
    _contours = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>(*source._contours));
    _holes = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>(*source._holes));
  }
  
  Polygon::Polygon(const SimplePolygon & simple_polygon)
  { 
    _contours = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>);
    _holes = std::auto_ptr< std::vector<SimplePolygon> >(new std::vector<SimplePolygon>);
  
    _contours->push_back(simple_polygon);
  }
  
  Polygon & Polygon::operator=(const Polygon & source)
  {
    *_contours = *source._contours;
    *_holes = *source._holes;
    return *this;
  }
  
  Polygon::Polygon(std::auto_ptr< std::vector<SimplePolygon> > contours, std::auto_ptr< std::vector<SimplePolygon> > & holes)
    : _contours(contours), _holes(holes) {}
  
  void Polygon::set_contours(std::auto_ptr< std::vector<SimplePolygon> > contours)
  {
    _contours.reset();
    _contours = contours;
  }
  
  void Polygon::set_holes(std::auto_ptr< std::vector<SimplePolygon> > holes)
  {
    _holes.reset();
    _holes = holes;
  }
 
 
  void Polygon::clear()
  {
    _contours->clear();
    _holes->clear();
  }
  
  
  /** \cond */
  void from_gpc_vertex_list(const gpc_vertex_list & vertex_list, SimplePolygon & contour)
  {
    std::auto_ptr<SimplePolygon::vertex_list_t> vertices(new SimplePolygon::vertex_list_t);
    
    vertices->resize(vertex_list.num_vertices);
    
    for(size_t j = 0; j < vertex_list.num_vertices; ++j)
    {
      (*vertices)[j](0) = vertex_list.vertex[j].x;
      (*vertices)[j](1) = vertex_list.vertex[j].y;
    }
    
    contour.set_vertices(vertices);
  }
  
  void to_gpc_polygon(const Polygon & polygon, gpc_polygon & out_gpc_poly)
  {
    gpc_vertex_list v_list;
    out_gpc_poly.num_contours = 0;
    out_gpc_poly.hole = 0;
    out_gpc_poly.contour = 0;
    
    for(size_t i = 0; i < polygon.n_contours(); ++i)
    {
      v_list.num_vertices = polygon.contour(i).n_vertices();
      v_list.vertex = (gpc_vertex*)&(polygon.contour(i).vertex(0));
      gpc_add_contour(&out_gpc_poly, &v_list, 0);
    }
    
    for(size_t i = 0; i < polygon.n_holes(); ++i)
    {
      v_list.num_vertices = polygon.hole(i).n_vertices();
      v_list.vertex = (gpc_vertex*)&(polygon.hole(i).vertex(0));
      gpc_add_contour(&out_gpc_poly, &v_list, 1);
    }
  }
  

  void from_gpc_poly(const gpc_polygon & in_gpc_poly, Polygon & polygon)
  {
    polygon.clear();
    
    size_t n_contours = 0;
    size_t n_holes = 0;
    
    std::auto_ptr< std::vector<SimplePolygon> > contours(new std::vector<SimplePolygon>);
    std::auto_ptr< std::vector<SimplePolygon> > holes(new std::vector<SimplePolygon>);
    
    for(size_t i = 0; i < in_gpc_poly.num_contours; ++i)
    {
      if(! in_gpc_poly.hole[i])
      {
        contours->resize(++n_contours);
        from_gpc_vertex_list(in_gpc_poly.contour[i], (*contours)[n_contours - 1]);
      }
      else
      {
        holes->resize(++n_holes);
        from_gpc_vertex_list(in_gpc_poly.contour[i], (*holes)[n_holes - 1]);
      }
    }
    
    polygon.set_contours(contours);
    polygon.set_holes(holes);
  }
  
  
  void polygon_union(const Polygon & poly_1, const Polygon & poly_2, Polygon & result)
  {
    gpc_polygon * gpc_poly_1 = new gpc_polygon;
    gpc_polygon * gpc_poly_2 = new gpc_polygon;
    
    to_gpc_polygon(poly_1, *gpc_poly_1);
    to_gpc_polygon(poly_2, *gpc_poly_2);
    
    gpc_polygon * gpc_result = new gpc_polygon;
    gpc_polygon_clip(GPC_UNION, gpc_poly_1, gpc_poly_2, gpc_result);
    
    from_gpc_poly(*gpc_result, result);
    
    gpc_free_polygon(gpc_poly_1);
    gpc_free_polygon(gpc_poly_2);
    gpc_free_polygon(gpc_result);
    
    delete gpc_poly_1;
    delete gpc_poly_2;
    delete gpc_result;
  }
  /** \endcond */
}

