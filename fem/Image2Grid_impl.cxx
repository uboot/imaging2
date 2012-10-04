#include <fem/Image2Grid_impl.hpp>

namespace imaging
{
  namespace image2grid_impl
  {
   ublas::fixed_vector<float_t, 2> mean_normal(const ublas::fixed_vector<float_t, 2> & reference, const ublas::fixed_vector<float_t, 2> & clockwise_neighbor, const ublas::fixed_vector<float_t, 2> & counter_clockwise_neighbor)
    {
      ublas::fixed_vector<float_t, 2> temp(clockwise_neighbor - counter_clockwise_neighbor);
      float_t norm = norm_2(temp);
      if(norm == 0.0)
        return ublas::fixed_vector<float_t, 2>(0.0, 0.0);
      else
        return ublas::fixed_vector<float_t, 2>(-temp(0), temp(1)) / norm_2(temp);
    }
    
    template <>
    void populate_grid<fem_2d_triangle_types>(Grid<fem_2d_triangle_types> & grid, const ublas::fixed_vector<size_t, 2> & size)
    {
      typedef Grid<fem_2d_triangle_types>::vertex_t vertex_t;
      
      size_t n_x_elements = size(0) - 1;
      size_t n_y_elements = size(1) - 1;
      
      for(size_t i = 0; i < n_x_elements; ++i)
        for(size_t j = 0; j < n_y_elements; ++j)
        {
          size_t element_index = i + j * n_x_elements;
          ublas::fixed_vector<size_t, 3> element_vertices_ll(element_index + j,
                                    element_index + j + 1,
                                    element_index + n_x_elements + j + 1);
          ublas::fixed_vector<size_t, 3> element_vertices_ur(element_index + n_x_elements + j + 2,
                                    element_index + n_x_elements + j + 1,
                                    element_index + j + 1);
                                    
          for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_vertices; ++m)
            grid.set_global_vertex_index(2 * element_index, m, element_vertices_ll(m));
          
          for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_nodes; ++m)
            grid.set_global_node_index(2 * element_index, m, element_vertices_ll(m));
          
          for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_vertices; ++m)
            grid.set_global_vertex_index(2 * element_index + 1, m, element_vertices_ur(m));
          
          for(size_t m = 0; m < Grid<fem_2d_triangle_types>::n_element_nodes; ++m)
            grid.set_global_node_index(2 * element_index + 1, m, element_vertices_ur(m));
  
  
          if( i == 0 )
          {
            grid.set_boundary_element(2 * (n_x_elements + n_y_elements) - i - j - 1,
                                      2 * element_index, 2);
          }
  
          if( i == n_x_elements - 1 )
          {
            grid.set_boundary_element(i + j + 1, 2 * element_index + 1, 2);
          }
  
          if( j == 0 )
          {
            grid.set_boundary_element(i + j, 2 * element_index, 0);
          }
  
          if( j == n_y_elements - 1 )
          {
            grid.set_boundary_element(2 * (n_x_elements + n_y_elements) - i - j - 2, 2 * element_index + 1, 0);
          }
        }
    }
    
    template <>
    void populate_grid<fem_2d_square_types>(Grid<fem_2d_square_types> & grid, const ublas::fixed_vector<size_t, 2> & size)
    {
      typedef Grid<fem_2d_square_types>::vertex_t vertex_t;
  
      size_t n_x_elements = size(0) - 1;
      size_t n_y_elements = size(1) - 1;
      
      for(size_t i = 0; i < n_x_elements; ++i)
        for(size_t j = 0; j < n_y_elements; ++j)
        {
          size_t element_index = i + j * n_x_elements;
          ublas::fixed_vector<size_t, 4> element_vertices(element_index + j,
                                    element_index + j + 1,
                                    element_index + n_x_elements + j + 2,
                                    element_index + n_x_elements + j + 1);
          
          
          for(size_t m = 0; m < Grid<fem_2d_square_types>::n_element_vertices; ++m)
            grid.set_global_vertex_index(element_index, m, element_vertices(m));
          
          for(size_t m = 0; m < Grid<fem_2d_square_types>::n_element_nodes; ++m)
            grid.set_global_node_index(element_index, m, element_vertices(m));
  
          if( i == 0 )
          {
            grid.set_boundary_element(2 * (n_x_elements + n_y_elements) - i - j - 1,
                                      element_index, 3);
          }
  
          if( i == n_x_elements - 1 )
          {
            grid.set_boundary_element(i + j + 1, element_index, 1);
          }
  
          if( j == 0 )
          {
            grid.set_boundary_element(i + j, element_index, 0);
          }
  
          if( j == n_y_elements - 1 )
          {
            grid.set_boundary_element(2 * (n_x_elements + n_y_elements) - i - j - 2, element_index, 2);
          }
        }
    }
    
    template <>
    void populate_grid<fem_3d_cube_types>(Grid<fem_3d_cube_types> & grid, const ublas::fixed_vector<size_t, 3> & size)
    {
      typedef Grid<fem_3d_cube_types>::vertex_t vertex_t;
      
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
      
      size_t n_x_elements = size(0) - 1;
      size_t n_y_elements = size(1) - 1;
      size_t n_z_elements = size(2) - 1;
      size_t boundary_index = 0;
      
      for(size_t i = 0; i < n_x_elements; ++i)
        for(size_t j = 0; j < n_y_elements; ++j)
          for(size_t k = 0; k < n_z_elements; ++k)
            
          {
            size_t element_index = i + j * n_x_elements + k * n_x_elements * n_y_elements;
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            
            ublas::fixed_vector<size_t, 8> element_vertices(vertex_index,
                                                            vertex_index + 1,
                                                            vertex_index + n_x_vertices,
                                                            vertex_index + n_x_vertices + 1,
                                                            vertex_index + n_x_vertices * n_y_vertices,
                                                            vertex_index + n_x_vertices * n_y_vertices + 1,
                                                            vertex_index + n_x_vertices * n_y_vertices + n_x_vertices,
                                                            vertex_index + n_x_vertices * n_y_vertices + n_x_vertices + 1);
            
            for(size_t m = 0; m < Grid<fem_3d_cube_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(element_index, m, element_vertices(m));
            
            for(size_t m = 0; m < Grid<fem_3d_cube_types>::n_element_nodes; ++m)
              grid.set_global_node_index(element_index, m, element_vertices(m));
            
            if( i == 0 )
            {
              grid.set_boundary_element(boundary_index, element_index, 2); 
              boundary_index = boundary_index + 1;
            }
            
            if( i == n_x_elements - 1 )
            {
              grid.set_boundary_element(boundary_index, element_index , 3); 
              boundary_index = boundary_index + 1;
            }
            
            if( j == 0 )
            {
              grid.set_boundary_element(boundary_index, element_index, 4); 
              boundary_index = boundary_index + 1;
            }
            
            if( j == n_y_elements - 1 )
            {
              grid.set_boundary_element(boundary_index, element_index, 5); 
              boundary_index = boundary_index + 1;
            }
            if( k == 0 )
            {
              grid.set_boundary_element(boundary_index, element_index, 0); 
              boundary_index = boundary_index + 1;
            }
            
            if( k == n_z_elements - 1 )
            {
              grid.set_boundary_element(boundary_index, element_index, 1); 
              boundary_index = boundary_index + 1;
            } 
            
          }
    }
    
    template <>
    void populate_grid<fem_3d_tetrahedra_types>(Grid<fem_3d_tetrahedra_types> & grid, const ublas::fixed_vector<size_t, 3> & size)
    {
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);    
      
      size_t n_x_elements = size(0) - 1;
      size_t n_y_elements = size(1) - 1;
      size_t n_z_elements = size(2) - 1;
      
      size_t boundary_index = 0;
      
      
      for(size_t i = 0; i < n_x_elements; ++i)
        for(size_t j = 0; j < n_y_elements; ++j)
          for(size_t k = 0; k < n_z_elements; ++k)
          {
            size_t element_index = i + j * n_x_elements + k * n_x_elements * n_y_elements;
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            
            ublas::fixed_vector<size_t, 4> element_vertices_one(vertex_index,
                                                                vertex_index + 1,
                                                                vertex_index + n_x_vertices,
                                                                vertex_index + n_x_vertices * n_y_vertices);
            ublas::fixed_vector<size_t, 4> element_vertices_two(vertex_index + 1,
                                                                vertex_index + n_x_vertices,
                                                                vertex_index + n_x_vertices * n_y_vertices,
                                                                vertex_index + n_x_vertices * n_y_vertices + n_x_vertices);
            ublas::fixed_vector<size_t, 4> element_vertices_three(vertex_index + 1,
                                                                  vertex_index + n_x_vertices * n_y_vertices,
                                                                  vertex_index + n_x_vertices * n_y_vertices +1,
                                                                  vertex_index + n_x_vertices * n_y_vertices + n_x_vertices);
            ublas::fixed_vector<size_t, 4> element_vertices_four(vertex_index + 1,
                                                                 vertex_index + n_x_vertices,
                                                                 vertex_index + n_x_vertices + 1,
                                                                 vertex_index + n_x_vertices * n_y_vertices + n_x_vertices);
            ublas::fixed_vector<size_t, 4> element_vertices_five(vertex_index + 1,
                                                                 vertex_index + n_x_vertices + 1,
                                                                 vertex_index + n_x_vertices * n_y_vertices + 1,
                                                                 vertex_index + n_x_vertices * n_y_vertices + n_x_vertices);
            ublas::fixed_vector<size_t, 4> element_vertices_six(vertex_index + n_x_vertices + 1,
                                                                vertex_index + n_x_vertices * n_y_vertices + 1,
                                                                vertex_index + n_x_vertices * n_y_vertices + n_x_vertices,
                                                                vertex_index + n_x_vertices * n_y_vertices + n_x_vertices +1);
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index, m, element_vertices_one(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index, m, element_vertices_one(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index + 1, m, element_vertices_two(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index + 1, m, element_vertices_two(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index + 2, m, element_vertices_three(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index + 2, m, element_vertices_three(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index + 3, m, element_vertices_four(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index + 3, m, element_vertices_four(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index + 4, m, element_vertices_five(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index + 4, m, element_vertices_five(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_vertices; ++m)
              grid.set_global_vertex_index(6 * element_index + 5, m, element_vertices_six(m));
            
            for(size_t m = 0; m < Grid<fem_3d_tetrahedra_types>::n_element_nodes; ++m)
              grid.set_global_node_index(6 * element_index + 5, m, element_vertices_six(m));
              
            if( i == 0 )
            {
              grid.set_boundary_element(boundary_index, 6 * element_index, 2); /* element_one */
              boundary_index = boundary_index + 1;
              
              grid.set_boundary_element(boundary_index, 6 * element_index + 1, 1); /* element_two */
              boundary_index = boundary_index + 1;
            }
            
            if( i == n_x_elements - 1 )
            {
              grid.set_boundary_element(boundary_index,6 * element_index + 4, 3); /* element_five */
              boundary_index = boundary_index + 1;
              grid.set_boundary_element(boundary_index,6 * element_index + 5, 0); /* element_six */
              boundary_index = boundary_index + 1;
            }
            
            if( j == 0 )
            {
              grid.set_boundary_element(boundary_index,6 * element_index, 0); /* element_one */
              boundary_index = boundary_index + 1;
              
              grid.set_boundary_element(boundary_index,6 * element_index + 2, 3); /* element_three */
              boundary_index = boundary_index + 1;
            }
            
            if( j == n_y_elements - 1 )
            {
              grid.set_boundary_element(boundary_index,6 * element_index + 3, 1); /* element_four */
              boundary_index = boundary_index + 1;
              
              grid.set_boundary_element(boundary_index,6 * element_index + 5, 2); /* element_six */
              boundary_index = boundary_index + 1;
            }
            if( k == 0 )
            {
              grid.set_boundary_element(boundary_index,6 * element_index, 3); /* element_one */
              boundary_index = boundary_index + 1;
              
              grid.set_boundary_element(boundary_index,6 * element_index + 3, 3); /* element_four */
              boundary_index = boundary_index + 1;
            }
            
            if( k == n_z_elements - 1 )
            {
              grid.set_boundary_element(boundary_index, 6 * element_index + 2, 1); /* element_three */
              boundary_index = boundary_index + 1;
              
              grid.set_boundary_element(boundary_index, 6 * element_index + 5, 1); /* element_six */
              boundary_index = boundary_index + 1;
            }
          }
    }
  }
}

