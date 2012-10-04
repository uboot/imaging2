#include <fem/utilities.hpp>

#include <core/utilities.hpp>
#include <shape/BoundaryDiscretizer.hpp>


extern "C"
{
  void triangle_triangulate(unsigned int n_in_points,
                   const double* in_points,
                   double max_triangle_area,
                   unsigned int *n_out_points,
                   double **out_points,
                   unsigned int *n_out_elements,
                   unsigned int **out_elements,
                   unsigned int *n_boundary_elements,
                   unsigned int **boundary_elements);
}


namespace imaging
{
  void uniform_grid(float_t lower_bound, float_t upper_bound, std::size_t n_elements, Grid<fem_1d_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size)
  {
    grid.set_dimensions(n_elements + 1, n_elements, 2, n_elements + 1);
    
    float_t element_length = (upper_bound - lower_bound) / n_elements;
    
    for(std::size_t i = 0; i <= n_elements; ++i)
      grid.set_vertex(i, Grid<fem_1d_types>::vertex_t(lower_bound + i * element_length));
      
    if(n_elements > 0)
    {
      grid.set_global_vertex_index(0, 0, 0);
      grid.set_global_vertex_index(0, 1, 1);
      grid.set_global_node_index(0, 0, 0);
      grid.set_global_node_index(0, 1, 1);
      //grid.set_boundary_node(0, ublas::fixed_vector<float_t, 1>(-1.0));
    }
    
    for(std::size_t i = 1; i < n_elements - 1; ++i)
    {
      grid.set_global_vertex_index(i, 0, i);
      grid.set_global_vertex_index(i, 1, i + 1);
      grid.set_global_node_index(i, 0, i);
      grid.set_global_node_index(i, 1, i + 1);
    }
    
    if(n_elements > 0)
    {
      grid.set_global_vertex_index(n_elements - 1, 0, n_elements - 1);
      grid.set_global_vertex_index(n_elements - 1, 1, n_elements);
      grid.set_global_node_index(n_elements - 1, 0, n_elements - 1);
      grid.set_global_node_index(n_elements - 1, 1, n_elements);
   //   grid.set_boundary_node(0, ublas::fixed_vector<float_t, 1>(1.0));
    }
    
    if(n_elements > 0)
    {
      grid.set_boundary_element(0, 0, 0);
      grid.set_boundary_element(1, n_elements - 1, 1);
    }
    stiffness_matrix_prototype.resize(system_size * (n_elements + 1), system_size * (n_elements + 1), false);
  }
  
  void circle_grid(float_t radius, std::size_t n_rings, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size)
  {
    ellipse_grid(radius, radius, n_rings, grid, stiffness_matrix_prototype, system_size);
  }
  
  
  void ellipse_grid(float_t a, float_t b, std::size_t n_rings, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size)
  {
    typedef ublas::fixed_vector<size_t, 3> element_vertices_t; 
    
    std::size_t n_vertices = 1 + 2 * (n_rings + 1) * n_rings;
    grid.set_dimensions(n_vertices, 
                        4 * n_rings * n_rings, 
                        4 * n_rings,
                        n_vertices);

    grid.set_vertex(0, Grid<fem_1d_types>::vertex_t(0.0, 0.0));
    std::size_t vertex = 1;
        
    for(std::size_t n = 1; n <= n_rings; ++n)
    {
      std::size_t current_n_vertices = 4 * n;
      float_t current_radius = float_t(n) / float_t(n_rings);
      float_t current_angle_offset = 2 * PI / float_t(current_n_vertices);
      
      for(std::size_t m = 0; m < current_n_vertices; ++m)
      {
        float_t current_angle = float_t(m) * current_angle_offset;
        grid.set_vertex(vertex, Grid<fem_1d_types>::vertex_t(a * current_radius * cos(current_angle),
                                                             b * current_radius * sin(current_angle)));
                                          
        if(n == n_rings)
        {
          ublas::fixed_vector<float_t, 2> unit_normal = ublas::fixed_vector<float_t, 2>(b * cos(current_angle), a * sin(current_angle));
          unit_normal /= norm_2(unit_normal);
          
         // grid.set_boundary_node(vertex, unit_normal);
        }
        
        vertex++;
      }
    }
    
    std::size_t start_index = 1;
    std::size_t element_index = 0;
    
    for(std::size_t n = 1; n <= n_rings; ++n)
    {
      std::size_t inner_start_index = start_index - max(4 * (n - 1), std::size_t(1));
      std::size_t outer_start_index = start_index + 4 * n;
      std::size_t inner_offset = 0;
      std::size_t outer_offset = 0;
      
      for(std::size_t m = 0; m < 4 * n; ++m)
      {
        ublas::fixed_vector<size_t, 3> element_vertices(start_index + m, inner_start_index + inner_offset % max(4 * (n - 1), std::size_t(1)), start_index + (m + 1) % (4 * n));
        
        for(size_t k = 0; k < Grid<fem_2d_triangle_types>::n_element_vertices; ++k)
          grid.set_global_vertex_index(element_index, k, element_vertices(k));
          
        for(size_t k = 0; k < Grid<fem_2d_triangle_types>::n_element_nodes; ++k)
          grid.set_global_node_index(element_index, k, element_vertices(k));
        
        if(n == n_rings)
        {
          grid.set_boundary_element(m, element_index, 2);
        }
      
        element_index++;
        
        if((4 * (m + 1)) % (4 * n) != 0)
         inner_offset++;
         
        if(n < n_rings)
        {
          element_vertices.assign(start_index + m, outer_start_index + (outer_offset + 1) % (4 * (n + 1)), start_index + (m + 1) % (4 * n));
          
          for(size_t k = 0; k < Grid<fem_2d_triangle_types>::n_element_vertices; ++k)
            grid.set_global_vertex_index(element_index, k, element_vertices(k));
          
          for(size_t k = 0; k < Grid<fem_2d_triangle_types>::n_element_nodes; ++k)
            grid.set_global_node_index(element_index, k, element_vertices(k));
          
          element_index++;
          
          outer_offset++;
          if((4 * (m + 1)) % (4 * n) == 0)
            outer_offset++;
        }
      }
      
      start_index += 4 * n;
    }
    stiffness_matrix_prototype.resize(system_size * n_vertices, system_size * n_vertices, false);
  }

  void triangulate_shape(const BoundaryDiscretizer<2> & shape_discretizer, float_t max_triangle_area, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size)
  {
    double *in_points = new double[2 * shape_discretizer.n_points()];
    unsigned int n_out_points;
    double *out_points;
    unsigned int n_out_elements;
    unsigned int *out_elements;
    unsigned int n_out_boundary_elements;
    unsigned int *out_boundary_elements;
    
    ublas::fixed_vector<float_t, 2> point, unit_normal;
    
    for(size_t i = 0; i < shape_discretizer.n_points(); ++i)
    {
      point = shape_discretizer(i, unit_normal);
      in_points[2 * i] = point(0);
      in_points[2 * i + 1] = point(1);
      
      unit_normal /= norm_2(unit_normal);
    //  grid.set_boundary_node(i, unit_normal);
    }
    
    triangle_triangulate(shape_discretizer.n_points(), in_points,
                         max_triangle_area,
                         &n_out_points, &out_points,
                         &n_out_elements, &out_elements,
                         &n_out_boundary_elements, &out_boundary_elements);
    
    grid.set_dimensions(n_out_points, n_out_elements, n_out_boundary_elements, n_out_points);
    
    for(size_t i = 0; i < grid.n_vertices(); ++i)
    {
      point(0) = out_points[2 * i];
      point(1) = out_points[2 * i + 1];
      
      grid.set_vertex(i, point);
    }
    
    ublas::fixed_vector<size_t, 3> element_vertices;
    
    for(size_t i = 0; i < grid.n_elements(); ++i)
    {
      element_vertices.assign(out_elements[3 * i], out_elements[3 * i + 1], out_elements[3 * i + 2]);
       
      for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_vertices; ++m)
        grid.set_global_vertex_index(i, m, element_vertices(m));
          
      for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_nodes; ++m)
        grid.set_global_node_index(i, m, element_vertices(m));
    }
    
    for(size_t i = 0; i < n_out_boundary_elements; ++i)
    {
      for(size_t j = 0; j < n_out_elements; ++j)
      {
        bool found_element = false;
        
        size_t k;
        for(k = 0; k < 3; ++k)
        {
          if( out_boundary_elements[2*i] == out_elements[3*j + k] && 
              out_boundary_elements[2*i + 1] == out_elements[3*j + (k + 1) % 3] )
            found_element = true;
          
          if(found_element)
            break;
          
          if( out_boundary_elements[2*i + 1] == out_elements[3*j + k] && 
              out_boundary_elements[2*i] == out_elements[3*j + (k + 1) % 3] )
            found_element = true;
          
          if(found_element)
            break;
        }
        
        if(found_element)
        {
          grid.set_boundary_element(i, j, k);
          break;
        }
      }
    }
    
    delete [] out_points;
    delete [] out_elements;
    delete [] out_boundary_elements;
    
    stiffness_matrix_prototype.resize(system_size * n_out_points, system_size * n_out_points, false);
  }
}

