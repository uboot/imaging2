// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_IMAGE2GRID_IMPL_H
#define FEM_IMAGE2GRID_IMPL_H

#include <fem/Grid.hpp>

#include <fem/fem_2d_square_types.hpp>
#include <fem/fem_2d_triangle_types.hpp>
#include <fem/fem_3d_cube_types.hpp>


namespace imaging
{
  namespace image2grid_impl
  {
    template <class fem_types>
    void populate_grid(Grid<fem_types> & grid, const ublas::fixed_vector<size_t, fem_types::data_dimension> & size);
    
    template <class vector_image_accessor_t, class fem_types>
    void construct_regular_2d_vertices(Grid<fem_types> & grid, const ublas::fixed_vector<size_t, 2> & size, const vector_image_accessor_t  & displacements)
    {
      typedef ublas::fixed_vector<float_t, 2> vertex_t;
      
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
  
      grid.set_regular(true);
  
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
        {
          size_t vertex_index = i + j * n_x_vertices;
          vertex_t offset(displacements[ublas::fixed_vector<size_t, 2>(i, j)]);
          grid.set_vertex(vertex_index, vertex_t( .5 + float_t(i) + offset(0), .5 + float_t(j) + offset(1)));
        }  
    }
    
    template <class vector_image_accessor_t, class fem_types>
    void construct_regular_3d_vertices(Grid<fem_types> & grid, const ublas::fixed_vector<size_t, 3> & size, const vector_image_accessor_t  & displacements)
    {
      typedef ublas::fixed_vector<float_t, 3> vertex_t;
      
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
      size_t n_z_vertices = size(2);
      
      grid.set_regular(false);
      
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for(size_t k = 0; k < n_z_vertices; ++k)
          {
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            vertex_t offset(displacements[ublas::fixed_vector<size_t, 3>(i, j, k)]);
            grid.set_vertex(vertex_index, vertex_t(0.5 + float_t(i) + offset(0), 
                                                   0.5 + float_t(j) + offset(1),
                                                   0.5 + float_t(k) + offset(2)));
          }  
    }
    
    template<class float_accessor_t> 
    void image2vector_2d(const ublas::fixed_vector<size_t, 2> size, const float_accessor_t & image, ublas::vector<float_t> & vector)
    {
      if(image.size() != size)
        throw Exception("Exception: Dimensions of image and grid do not agree in Image2Grid::image2vector().");

      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
  
      vector.resize(n_x_vertices * n_y_vertices, false);
  
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
        {
          size_t vertex_index = i + j * n_x_vertices;
          vector(vertex_index) = image[ublas::fixed_vector<size_t, 2>(i, j)];
        }
    }
     
    template<class float_accessor_t> 
    void image2boundary_vector_2d(const ublas::fixed_vector<size_t, 2> size, const float_accessor_t &image, ublas::mapped_vector<float_t> & vector)
    {
      if(image.size() != size)
        throw Exception("Exception: Dimensions of image and grid do not agree in Image2Grid::image2boundary_vector().");
  
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
  
      vector.resize(n_x_vertices * n_y_vertices, false);
  
      for(size_t i = 0; i < n_x_vertices; ++i)
        for(size_t j = 0; j < n_y_vertices; ++j)
        {
          size_t vertex_index = i + j * n_x_vertices;
          vector(vertex_index) = image[ublas::fixed_vector<size_t, 2>(i, j)];
  
          if(i > 0 && i < n_x_vertices - 1)
            j += n_y_vertices - 2;
        }
    }
    
    template <class float_accessor_t>
    void vector2image_2d(const ublas::fixed_vector<size_t, 2> size, const ublas::vector< float_t > &vector, float_accessor_t & image)
    {
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
  
      if(vector.size() != n_x_vertices * n_y_vertices)
        throw Exception("Exception: vector of wrong length in vector2image().");
      
      if(image.size() != size)
        throw Exception("Exception: image of wrong size in vector2image().");
  
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
        {
          size_t vertex_index = i + j * n_x_vertices;
          image[ublas::fixed_vector<size_t, 2>(i, j)] = vector(vertex_index);
        }
    }
    
    template <class vector_image_accessor_t>
    void vector2vector_image_2d(const ublas::fixed_vector<size_t, 2> size, const ublas::vector< float_t > &vector, vector_image_accessor_t & vector_image)
    {
      const size_t n_x_vertices = size(0);
      const size_t n_y_vertices = size(1);
      const size_t dimension = vector_image_accessor_t::data_t::dimension;
  
      if(vector.size() != n_x_vertices * n_y_vertices * dimension)
        throw Exception("Exception: tla::vector of wrong length in vector2image().");
      
      if(vector_image.size() != size)
        throw Exception("Exception: vector_image of wrong size in vector2vector_image().");
  
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for(size_t k = 0; k < dimension; ++k)
          {
            size_t vertex_index = dimension * ( i + j * n_x_vertices ) + k;
            vector_image[ublas::fixed_vector<size_t, 2>(i, j)](k) = vector(vertex_index);
          }
    }
    
    template<class float_accessor_t> 
    void image2vector_3d(const ublas::fixed_vector<size_t, 3> size, const float_accessor_t & image, ublas::vector<float_t> & vector)
    {
      if(image.size() != size)
        throw Exception("Exception: Dimensions of image and grid do not agree in Image2Grid::image2vector().");
      
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
      size_t n_z_vertices = size(2);
      
      vector.resize(n_x_vertices * n_y_vertices * n_z_vertices, false);
      
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for (size_t k = 0; k < n_z_vertices; ++k)
          {
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            vector(vertex_index) = image[ublas::fixed_vector<size_t, 3>(i, j, k)];
          }
    }
     
    template<class float_accessor_t> 
    void image2boundary_vector_3d(const ublas::fixed_vector<size_t, 3> size, const float_accessor_t &image, ublas::mapped_vector<float_t> & vector)
    {
      if(image.size() != size)
        throw Exception("Exception: Dimensions of image and grid do not agree in Image2Grid::image2boundary_map().");
      
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
      size_t n_z_vertices = size(2);
     
      vector.resize(n_x_vertices * n_y_vertices * n_z_vertices, false);
      
      for(size_t i = 0; i < n_x_vertices; ++i)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for(size_t k = 0; k < n_z_vertices; ++k)  
          {
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            vector(vertex_index) = image[ublas::fixed_vector<size_t, 3>(i, j, k)];
            if(i > 0 && i < n_x_vertices - 1 && j > 0 && j < n_y_vertices -1 )
              k += n_z_vertices - 2;
          }
    }
    
    template <class float_accessor_t>
    void vector2image_3d(const ublas::fixed_vector<size_t, 3> size, const ublas::vector< float_t > &vector, float_accessor_t & image)
    {
      size_t n_x_vertices = size(0);
      size_t n_y_vertices = size(1);
      size_t n_z_vertices = size(2);
      
      if(vector.size() != n_x_vertices * n_y_vertices * n_z_vertices)
        throw Exception("Exception: tla::vector of wrong length in vector2image().");
      
      if(image.size() != size)
        throw Exception("Exception: image of wrong size in vector2image().");
      
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for(size_t k = 0; k < n_z_vertices; ++k) 
          {
            size_t vertex_index = i + j * n_x_vertices + k * n_x_vertices * n_y_vertices;
            image[ublas::fixed_vector<size_t, 3>(i, j, k)] = vector(vertex_index);
          }
    }
    
    template <class vector_image_accessor_t>
    void vector2vector_image_3d(const ublas::fixed_vector<size_t, 3> size, const ublas::vector< float_t > &vector, vector_image_accessor_t & vector_image)
    {
      const size_t n_x_vertices = size(0);
      const size_t n_y_vertices = size(1); 
      const size_t n_z_vertices = size(2);
      const size_t dimension = vector_image_accessor_t::data_t::dimension;
      
      if(vector.size() != n_x_vertices * n_y_vertices * n_z_vertices * dimension)
        throw Exception("Exception: tla::vector of wrong length in vector2image().");
      
      if(vector_image.size() != size)
        throw Exception("Exception: vector_image of wrong size in vector2vector_image().");
      
      for(size_t i = 0; i < n_x_vertices; i++)
        for(size_t j = 0; j < n_y_vertices; ++j)
          for(size_t k = 0; k < n_z_vertices; ++k)  
            for(size_t l = 0; l < dimension; ++l)
            {
              size_t vertex_index = dimension * ( i + j * n_x_vertices + k * n_x_vertices * n_y_vertices) + l;
              vector_image[ublas::fixed_vector<size_t, 3>(i, j, k)](l) = vector(vertex_index);
            }
    }
    
    
    template <>
    void populate_grid<fem_2d_triangle_types>(Grid<fem_2d_triangle_types> & grid, const ublas::fixed_vector<size_t, 2> & size);
    
    template <>
    void populate_grid<fem_2d_square_types>(Grid<fem_2d_square_types> & grid, const ublas::fixed_vector<size_t, 2> & size);

	  template <>
    void populate_grid<fem_3d_cube_types>(Grid<fem_3d_cube_types> & grid, const ublas::fixed_vector<size_t, 3> & size);
    
    template <>
    void populate_grid<fem_3d_tetrahedra_types>(Grid<fem_3d_tetrahedra_types> & grid, const ublas::fixed_vector<size_t, 3> & size);
  }
}


#endif
