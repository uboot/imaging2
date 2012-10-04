#include <fem/Image2Grid.hpp>


namespace imaging
{
  /** \cond */
  template <>
  void Image2Grid<fem_2d_square_types>::stiffness_matrix_prototype(ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, size_t system_size) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);
    size_t n_nodes = n_x_vertices * n_y_vertices;
    
    ublas::compressed_matrix<float_t> & matrix = stiffness_matrix_prototype;
    matrix.resize(system_size * n_nodes, system_size * n_nodes, false);

    for(size_t j = 0; j < n_nodes; ++j)
    {
      for(size_t i = system_size * j; i < system_size * (j + 1); ++i)
      {
        if(j >= n_x_vertices + 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * ( n_x_vertices + 1 ), 0.0);
        if(j >= n_x_vertices) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * n_x_vertices, 0.0);
        if(j >= n_x_vertices - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * ( n_x_vertices - 1 ), 0.0);
  
        if(j >= 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size, 0.0);
        for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
          matrix.insert_element(i, k, 0.0);
        if(j < n_nodes - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size, 0.0);
  
        if(j < n_nodes - n_x_vertices + 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * ( n_x_vertices - 1 ), 0.0);
        if(j < n_nodes - n_x_vertices) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * n_x_vertices, 0.0);
        if(j < n_nodes - n_x_vertices - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * ( n_x_vertices + 1 ), 0.0);
      }
    }
  }

  template <>
  void Image2Grid<fem_2d_triangle_types>::stiffness_matrix_prototype(ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, size_t system_size) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);
    size_t n_nodes = n_x_vertices * n_y_vertices;
    
    ublas::compressed_matrix<float_t> & matrix = stiffness_matrix_prototype;
    matrix.resize(system_size * n_nodes, system_size * n_nodes, false);

    for(size_t j = 0; j < n_nodes; ++j)
    {
      for(size_t i = system_size * j; i < system_size * (j + 1); ++i)
      {
        if(j >= n_x_vertices + 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * ( n_x_vertices + 1 ), 0.0);
        if(j >= n_x_vertices) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * n_x_vertices, 0.0);
        if(j >= n_x_vertices - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size * ( n_x_vertices - 1 ), 0.0);
  
        if(j >= 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k - system_size, 0.0);
        for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
          matrix.insert_element(i, k, 0.0);
        if(j < n_nodes - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size, 0.0);
  
        if(j < n_nodes - n_x_vertices + 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * ( n_x_vertices - 1 ), 0.0);
        if(j < n_nodes - n_x_vertices) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * n_x_vertices, 0.0);
        if(j < n_nodes - n_x_vertices - 1) 
          for(size_t k = system_size * j; k < system_size * (j + 1); ++k)
            matrix.insert_element(i, k + system_size * ( n_x_vertices + 1 ), 0.0);
      }
    }
    
    
  }
  
  /** \endcond */

}

