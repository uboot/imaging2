#include <polytope/SimplePolygon.hpp>

extern "C" {
double inter(const double * a, int na, const double * b, int nb);
}

namespace imaging
{
  SimplePolygon::SimplePolygon()
  {
    _vertices = std::auto_ptr<vertex_list_t>(new vertex_list_t);
  }
  
  SimplePolygon::SimplePolygon(const SimplePolygon & source)
  {
    _vertices = std::auto_ptr<vertex_list_t>(new vertex_list_t(*source._vertices));
  }
  
  SimplePolygon & SimplePolygon::operator=(const SimplePolygon & source)
  {
    *_vertices = *source._vertices;
    return *this;
  }
  
  /** Returns the <em>i</em>-th vertex of the polygon. */
  const ublas::fixed_vector<float_t, 2> & SimplePolygon::vertex(size_t i) const
  { 
    return (*_vertices)[i];
  }
  
  /** Returns the total number of vertices of the polygon. */
  size_t SimplePolygon::n_vertices() const 
  {
    return _vertices->size();
  }
  
  void SimplePolygon::set_vertices(std::auto_ptr<vertex_list_t> vertices)
  {
    _vertices.reset();
    _vertices = vertices;
  }
  
  /** Returns the vertices of the polygon. */
  const SimplePolygon::vertex_list_t & SimplePolygon::vertices() const { return *_vertices; }
    
  float_t compute_intersection_volume(const SimplePolygon & poly_a, const SimplePolygon & poly_b)
  {
    return inter(&(poly_a.vertex(0)(0)), poly_a.n_vertices(), &(poly_b.vertex(0)(0)), poly_b.n_vertices());
  }
}
