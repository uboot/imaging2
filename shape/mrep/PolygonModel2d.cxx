#include <shape/mrep/PolygonModel2d.hpp>

#include <core/utilities.hpp>
#include <core/matrix_utilities.hpp>
#include <core/vector_utilities.hpp>

namespace imaging
{
  void discretize_sector(ublas::fixed_vector<float_t, 2> center, float_t radius, float_t angle_1, float_t ccw_offset, size_t n_vertices,
                         std::vector< ublas::fixed_vector<float_t, 2> > vertices)
  {
    vertices.resize(n_vertices);
    float_t offset = ccw_offset / float_t(n_vertices);
    
    for(size_t i = 0; i < n_vertices; ++i)
      vertices[i] = center + polar2cartesian(radius, angle_1 + i * offset);
  }
  
  
  void PolygonModel2d::compute_limb_polygon(size_t atom_1, size_t atom_2, size_t n_sector_discretization_points, 
    Polygon & polygon) const
  {
    n_sector_discretization_points++;
    
    ublas::fixed_vector<float_t, 2> atom_center_1(atom_center(atom_1));
    ublas::fixed_vector<float_t, 2> atom_center_2(atom_center(atom_2));
     
    float_t normal_angle_1 = compute_tangent_angle(atom_2, atom_1);
    float_t normal_angle_2 = compute_tangent_angle(atom_1, atom_2);
    
    float_t sector_angle_2 = counter_clockwise_difference(normal_angle_1, normal_angle_2);
    float_t sector_angle_1 = 2 * PI - sector_angle_2;
    
    ublas::fixed_vector<float_t, 2> normal_1(polar2cartesian (1.0, normal_angle_1));
    ublas::fixed_vector<float_t, 2> normal_2(polar2cartesian (1.0, normal_angle_2));
    
    std::auto_ptr< std::vector<SimplePolygon> > contours(new std::vector<SimplePolygon>(1));
    std::auto_ptr<SimplePolygon::vertex_list_t> vertices(new SimplePolygon::vertex_list_t);

    vertices->push_back(atom_center_1 + atom_radius(atom_1) * normal_1);
    
    float_t offset = sector_angle_2 / float_t(n_sector_discretization_points);
    
    for(size_t i = 0; i < n_sector_discretization_points; ++i)
      vertices->push_back(atom_center_2 + polar2cartesian(atom_radius(atom_2), normal_angle_1 + float_t(i) * offset));
    
    vertices->push_back(atom_center_2 + atom_radius(atom_2) * normal_2);
   
    offset = sector_angle_1 / float_t(n_sector_discretization_points);
    
    for(size_t i = 0; i < n_sector_discretization_points; ++i)
      vertices->push_back(atom_center_1 + polar2cartesian(atom_radius(atom_1), normal_angle_2 + float_t(i) * offset));
    
    polygon.clear();
    (*contours)[0].set_vertices(vertices);
    polygon.set_contours(contours);
  }
  
  
  void PolygonModel2d::compute_boundary(size_t n_sector_discretization_points, Polygon & polygon) const
  {
    polygon.clear();
    Polygon limb_polygon;
      
    for(size_t i = 1; i < n_atoms(); ++i)
    {
      compute_limb_polygon(atom(i).start_atom(), i, n_sector_discretization_points, limb_polygon);
      polygon_union(polygon, limb_polygon, polygon);
    }
  }
  
  
  PolygonModel2d::Discretizer::Discretizer(const PolygonModel2d & model, size_t n_points) : img::BoundaryDiscretizer<2>(n_points)
  {
    Polygon boundary_polygon;
    model.compute_boundary(N_SECTOR_DISCRETIZATION_POINTS, boundary_polygon);
    
    float_t boundary_length = 0.0;
    
    for(size_t i = 0; i < boundary_polygon.contours().size(); ++i)
      for(size_t j = 0; j < boundary_polygon.contours()[i].n_vertices(); ++j)
        boundary_length += norm_2(boundary_polygon.contours()[i].vertex( (j + 1) % boundary_polygon.contours()[i].n_vertices()) - 
                                  boundary_polygon.contours()[i].vertex(j));
    
    for(size_t i = 0; i < boundary_polygon.holes().size(); ++i)
      for(size_t j = 0; j < boundary_polygon.holes()[i].n_vertices(); ++j)
        boundary_length += norm_2(boundary_polygon.holes()[i].vertex( (j + 1) % boundary_polygon.holes()[i].n_vertices()) - 
                                  boundary_polygon.holes()[i].vertex(j));
    
    float_t step_size = boundary_length / n_points;
    float_t current_length = step_size;
    
    for(size_t i = 0; i < boundary_polygon.contours().size(); ++i)
      for(size_t j = 0; j < boundary_polygon.contours()[i].n_vertices(); ++j)
      {
        const ublas::fixed_vector<float_t, 2> & first_point = boundary_polygon.contours()[i].vertex(j);
        const ublas::fixed_vector<float_t, 2> & second_point = boundary_polygon.contours()[i].vertex( (j + 1) % boundary_polygon.contours()[i].n_vertices());
        
        ublas::fixed_vector<float_t, 2> segment_vector(second_point - first_point);
        float_t segment_length = norm_2(segment_vector);
        float_t remaining_segment = segment_length;
        float_t used_segment = 0.0;
        
        while(remaining_segment > current_length)
        {
          _points.push_back(first_point + ( current_length + used_segment ) / segment_length * segment_vector);
          _normals.push_back(step_size / norm_2(segment_vector) * 
                             ublas::fixed_vector<float_t, 2>(- segment_vector(1), segment_vector(0)));
          used_segment += current_length;
          remaining_segment -= current_length;
          current_length = step_size;
        }
        
        current_length -= remaining_segment;
      }
    
    for(size_t i = 0; i < boundary_polygon.holes().size(); ++i)
      for(size_t j = 0; j < boundary_polygon.holes()[i].n_vertices(); ++j)
      {
        const ublas::fixed_vector<float_t, 2> & first_point = boundary_polygon.holes()[i].vertex(j);
        const ublas::fixed_vector<float_t, 2> & second_point = boundary_polygon.holes()[i].vertex( (j + 1) % boundary_polygon.holes()[i].n_vertices());
        
        ublas::fixed_vector<float_t, 2> segment_vector(second_point - first_point);
        float_t segment_length = norm_2(segment_vector);
        float_t remaining_segment = segment_length;
        float_t used_segment = 0.0;
        
        while(remaining_segment > current_length)
        {
          _points.push_back(first_point + ( current_length + used_segment ) / segment_length * segment_vector);
          _normals.push_back(step_size / norm_2(segment_vector) * 
                             ublas::fixed_vector<float_t, 2>(- segment_vector(1), segment_vector(0)));
          used_segment += current_length;
          remaining_segment -= current_length;
          current_length = step_size;
        }
        
        current_length -= remaining_segment;
      }
        
  }

  void PolygonModel2d::Discretizer::evaluate(size_t i, ublas::fixed_vector<float_t, 2> & point, ublas::fixed_vector<float_t, 2> & normal, float_t & curvature) const
  {
    if(i < _points.size())
    {
      point = _points[i];
      normal = _normals[i];
    }
  }
  
  std::auto_ptr< BoundaryDiscretizer<2> > PolygonModel2d::boundary_discretizer(size_t n_points) const
  {
    return std::auto_ptr< BoundaryDiscretizer<2> >(new Discretizer(*this, n_points));
  }

}

