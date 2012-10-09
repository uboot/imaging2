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

#ifndef FEM_IMAGE2GRID_H
#define FEM_IMAGE2GRID_H

#include <fem/Image2Grid_impl.hpp>
#include <image/Image.hpp>
#include <image/ScalarImage.hpp>


namespace imaging
{
  /** \ingroup fem
      \brief Constructs FE grids from multi-dimensional images and provides conversion functions between grid and image.
      
       This class constructs a FE grid for a given image which reflects the geomotry of the image data. The type of the grid is determined by the template parameter \em fem_types. Its dimension (\em fem_types::data_dimension) must be the same as the dimension of input image data. Additionally the class provides functions to convert (multi-dimensional) image data to the 1-dimensional vector representation necessary for FE computations and back again.
       
       Note that all the functions in this class are implemented for specific \em fem_types classes only (mainly fem_2d_square_types and fem_2d_triangle_types).
       For \em fem_2d_square_types the structure of a pixel is the following:
       \f[
         \begin{array}{ccc}
           3 & - & 2 \\
           | &   & | \\
           0 & - & 1
         \end{array}
       \f]
       For \em fem_2d_triangle_types the structure of a pixel is the following:
       \f[
         \begin{array}{ccc}
           2 &   &  \\
           | &   &  \\
           0 & - & 1
         \end{array}
         \quad \textrm{and} \quad
         \begin{array}{ccc}
           1 & - & 0 \\
             &   & | \\
             &   & 2
         \end{array}
       \f]
       
  */
  template<class fem_types>
  class Image2Grid
  {
    typedef ublas::fixed_vector<float_t, fem_types::data_dimension> vertex_t;
    ublas::fixed_vector<size_t, fem_types::data_dimension> _size;
    
  public:
    /** Constructs an Image2Grid objects which creates the grid determined by \em dimensions or converts data from and to it. */
    Image2Grid(const ublas::fixed_vector<size_t, fem_types::data_dimension> & dimensions) { _size = dimensions; }

    /** Sets the geometry of \em grid to the dimensions specified in the constructor of the Image2Grid object. */
    void construct_grid(Grid<fem_types> & grid) const  
    {
      ublas::fixed_vector<float_t, fem_types::data_dimension> zeros(0.0);
      construct_grid(grid, ScalarImage<fem_types::data_dimension, ublas::fixed_vector<float_t, fem_types::data_dimension> >(_size, zeros));
    }

    /** Sets the geometry of \em grid to the dimensions specified in the constructor of the Image2Grid object. Each node of the grid is displaced by the displacement vectors in \em displacements. */
    template <class vector_image_accessor_t>
    void construct_grid(Grid<fem_types> & grid, const vector_image_accessor_t & displacements) const;

    /** Resizes \em stiffness_matrix to the correct size for \em grid and prefilled with values at matrix position where non-zero entries are expected (this depends on the geometry of the grid). In case the PDE to be solved is not scalar but a system of equation, the user has to pass the number of equations (\em system_size) to ensure that \em stiffness_matrix is sized correctly. */
    void stiffness_matrix_prototype(ublas::compressed_matrix<float_t> & stiffness_matrix, size_t system_size = 1) const
    {
      size_t n_nodes = 1;
      
      for(size_t i = 0; i < fem_types::data_dimension; ++i)
        n_nodes *= _size(i);
    
      stiffness_matrix.resize(n_nodes * system_size, n_nodes * system_size, false);
    }

    /** Converts \em image to \em vector such that the indices of the values in \em vector conform with the node indices of the grid constructed by this Image2Grid object.
        \param[in] image An Image or ImageAccessorInterface object which has scalar (grayscale) values. The dimensions of \em image must match the dimensions of the grid, otherwise an Exception is thrown.
        \param[out] vector Is automatically resized to the number of nodes in the grid.
    */
    template <class float_accessor_t>
    void image2vector(const float_accessor_t & image, ublas::vector<float_t> & vector) const;

    /** Converts \em vector to \em image such that the indices of the values in \em vector conform with the node indices of the grid constructed by this Image2Grid object.
        \param[in] vector Its length must match the number of nodes on the grid, otherwise an Exception is thrown. 
        \param[out] image An Image or ImageAccessorInterface object which has scalar (grayscale) values. The dimensions of \em image must match the dimensions of the grid, otherwise an Exception is thrown.
    */
    template <class float_accessor_t>
    void vector2image(const ublas::vector<float_t> & vector, float_accessor_t & image) const;
    
    /** Converts \em vector to \em vector_image such that the indices of the values in \em vector conform with the node indices of the grid constructed by this Image2Grid object.
        \param[in] vector Its length must match the number of nodes on the grid times the dimension of the values of \em vector_image, otherwise an Exception is thrown. 
        \param[out] vector_image An Image or ImageAccessorInterface object which has vector values. The dimensions of \em vector_image must match the dimensions of the grid, otherwise an Exception is thrown.
    */    
    template <class vector_image_accessor_t>
    void vector2vector_image(const ublas::vector<float_t> & vector, vector_image_accessor_t & vector_image) const;

    /** Converts \em image to \em vector such that the indices of the values in \em vector conform with the node indices of the grid constructed by this Image2Grid object. In contrast to image2vector() only the image values at the boundary are copied.
        \param[in] image An Image or ImageAccessorInterface object which has scalar (grayscale) values. The dimensions of \em image must match the dimensions of the grid, otherwise an Exception is thrown.
        \param[out] vector Is automatically resized to the number of nodes in the grid.
    */
    template <class float_accessor_t>
    void image2boundary_vector(const float_accessor_t & image, ublas::mapped_vector<float_t> & vector)const;
  }
  ;
  
  /** \cond */
  template <>
  void Image2Grid<fem_2d_square_types>::stiffness_matrix_prototype(ublas::compressed_matrix<float_t> & stiffness_matrix, size_t system_size) const;

  template <>
  void Image2Grid<fem_2d_triangle_types>::stiffness_matrix_prototype(ublas::compressed_matrix<float_t> & stiffness_matrix, size_t system_size) const;
  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_square_types>::image2vector(const float_accessor_t & image, ublas::vector<float_t> & vector) const
  {
    image2grid_impl::image2vector_2d(_size, image, vector);
  }


  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_square_types>::vector2image(const ublas::vector<float_t> & vector, float_accessor_t & image) const
  {
    image2grid_impl::vector2image_2d(_size, vector, image);
  }


  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_2d_square_types>::vector2vector_image(const ublas::vector<float_t> & vector, vector_image_accessor_t & vector_image) const
  {
    image2grid_impl::vector2vector_image_2d(_size, vector, vector_image);
  }


  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_square_types>::image2boundary_vector(const float_accessor_t & image, ublas::mapped_vector<float_t> & vector) const
  {
    image2grid_impl::image2boundary_vector_2d(_size, image, vector);
  }
  
    template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_triangle_types>::image2vector(const float_accessor_t & image, ublas::vector<float_t> & vector) const
  {
    image2grid_impl::image2vector_2d(_size, image, vector);
  }


  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_triangle_types>::vector2image(const ublas::vector<float_t> & vector, float_accessor_t & image) const
  {
    image2grid_impl::vector2image_2d(_size, vector, image);
  } 
  
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_2d_triangle_types>::vector2vector_image(const ublas::vector<float_t> & vector, vector_image_accessor_t & vector_image) const
  {
    image2grid_impl::vector2vector_image_2d(_size, vector, vector_image);
  }


  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_2d_triangle_types>::image2boundary_vector(const float_accessor_t & image, ublas::mapped_vector<float_t> & vector) const
  {
    image2grid_impl::image2boundary_vector_2d(_size, image, vector);
  }


  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_cube_types>::image2vector(const float_accessor_t & image, ublas::vector<float_t> & vector) const
  {
    image2grid_impl::image2vector_3d(_size, image, vector);
  }
  
  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_cube_types>::vector2image(const ublas::vector<float_t> & vector, float_accessor_t & image) const
  {
    image2grid_impl::vector2image_3d(_size, vector, image);
  } 
  
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_3d_cube_types>::vector2vector_image(const ublas::vector<float_t> & vector, vector_image_accessor_t & vector_image) const
  {
    image2grid_impl::vector2vector_image_3d(_size, vector, vector_image);
  }
  
  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_cube_types>::image2boundary_vector(const float_accessor_t & image, ublas::mapped_vector<float_t> & vector) const
  {
    image2grid_impl::image2boundary_vector_3d(_size, image, vector);
  }

  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_tetrahedra_types>::image2vector(const float_accessor_t & image, ublas::vector<float_t> & vector) const
  {
    image2grid_impl::image2vector_3d(_size, image, vector);
  }
  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_tetrahedra_types>::vector2image(const ublas::vector<float_t> & vector, float_accessor_t & image) const
  {
    image2grid_impl::vector2image_3d(_size, vector, image);
  } 
  
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_3d_tetrahedra_types>::vector2vector_image(const ublas::vector<float_t> & vector, vector_image_accessor_t & vector_image) const
  {
    image2grid_impl::vector2vector_image_3d(_size, vector, vector_image);
  }
  
  
  template<>
  template <class float_accessor_t>
  void Image2Grid<fem_3d_tetrahedra_types>::image2boundary_vector(const float_accessor_t & image, ublas::mapped_vector<float_t> & vector) const
  {
    image2grid_impl::image2boundary_vector_3d(_size, image, vector);
  }
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_2d_square_types>::construct_grid(Grid<fem_2d_square_types> & grid, const vector_image_accessor_t  & displacements) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);
    size_t n_x_elements = _size(0) - 1;
    size_t n_y_elements = _size(1) - 1;

    grid.set_dimensions( n_x_vertices * n_y_vertices,
                         n_x_elements * n_y_elements,
                         2 * (n_x_elements + n_y_elements),
                         n_x_vertices * n_y_vertices );
                         
    grid.set_regular(true);
                         
    image2grid_impl::construct_regular_2d_vertices(grid, _size, displacements);
    image2grid_impl::populate_grid(grid, _size);
  }
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_2d_triangle_types>::construct_grid(Grid<fem_2d_triangle_types> & grid, const vector_image_accessor_t  & displacements) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);
    size_t n_x_elements = _size(0) - 1;
    size_t n_y_elements = _size(1) - 1;

    grid.set_dimensions( n_x_vertices * n_y_vertices,
                         2 * (n_x_elements * n_y_elements),
                         2 * (n_x_elements + n_y_elements),
                         n_x_vertices * n_y_vertices );
                         
    grid.set_regular(true);
                         
    image2grid_impl::construct_regular_2d_vertices(grid, _size, displacements);
    image2grid_impl::populate_grid(grid, _size);
  }
  
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_3d_cube_types>::construct_grid(Grid<fem_3d_cube_types> & grid, const vector_image_accessor_t  & displacements) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);    
    size_t n_z_vertices = _size(2);
    size_t n_x_elements = _size(0) - 1;
    size_t n_y_elements = _size(1) - 1;
    size_t n_z_elements = _size(2) - 1;
    
    grid.set_dimensions(n_x_vertices * n_y_vertices * n_z_vertices, /* number of vertices*/
                        n_x_elements * n_y_elements * n_z_elements,
                        2 * (n_x_elements * n_y_elements + n_y_elements * n_z_elements + n_z_elements * n_x_elements),
                        n_x_vertices * n_y_vertices * n_z_vertices);
                        
    grid.set_regular(true);
    
    image2grid_impl::construct_regular_3d_vertices(grid, _size, displacements);
    image2grid_impl::populate_grid(grid, _size);
  }
  
    
  template<>
  template <class vector_image_accessor_t>
  void Image2Grid<fem_3d_tetrahedra_types>::construct_grid(Grid<fem_3d_tetrahedra_types> & grid, const vector_image_accessor_t  & displacements) const
  {
    size_t n_x_vertices = _size(0);
    size_t n_y_vertices = _size(1);    
    size_t n_z_vertices = _size(2);
    size_t n_x_elements = _size(0) - 1;
    size_t n_y_elements = _size(1) - 1;
    size_t n_z_elements = _size(2) - 1;
    
    grid.set_dimensions( n_x_vertices * n_y_vertices * n_z_vertices, /* number of vertices*/
                        6 * (n_x_elements * n_y_elements * n_z_elements),
                        4 * (n_x_elements * n_y_elements + n_y_elements * n_z_elements + n_z_elements * n_x_elements),
                        n_x_vertices * n_y_vertices * n_z_vertices);
                        
    grid.set_regular(false);
    
    image2grid_impl::construct_regular_3d_vertices(grid, _size, displacements);
    image2grid_impl::populate_grid(grid, _size);
  }
  /** \endcond */
}


#endif
