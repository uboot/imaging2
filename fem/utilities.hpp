// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_UTILITIES_H
#define FEM_UTILITIES_H

#include <fem/Grid.hpp>
#include <fem/fem_2d_triangle_types.hpp>
#include <fem/fem_1d_types.hpp>

namespace imaging
{
  /** \ingroup fem
      <tt>\#include <fem/utilities.hpp></tt>
      
      Constructs a regular \em grid (i.e. with elements of constant size) on the interval defined by \em lower_bound and \em upper_bound for a system of equations. In addition, \em stiffness_matrix_prototype is resized to the correct size for \em grid. In case the PDE to be solved is not scalar but a system of equation, the user has to pass the number of equations (\em system_size) to ensure that \em stiffness_matrix_prototype is sized correctly.  
  */
  void uniform_grid(float_t lower_bound, float_t upper_bound, std::size_t n_elements, Grid<fem_1d_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size = 1);
  
  /** \ingroup fem
      <tt>\#include <fem/utilities.hpp></tt>
      
      Constructs a circular \em grid of \em radius. The triangular elements are laid out in \em n_rings rings. The <em>i</em>-th ring contains 4 * \em i elements. 
      In addition, \em stiffness_matrix_prototype is resized to the correct size for \em grid. In case the PDE to be solved is not scalar but a system of equation, the user has to pass the number of equations (\em system_size) to ensure that \em stiffness_matrix_prototype is sized correctly.  
  */
  void circle_grid(float_t radius, std::size_t n_rings, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size = 1);
  
    /** \ingroup fem
      <tt>\#include <fem/utilities.hpp></tt>
      
      Constructs a elliptic \em grid with axis of length \em a and \em b. The triangular elements are laid out in \em n_rings rings. The <em>i</em>-th ring contains 4 * \em i elements. 
      In addition, \em stiffness_matrix_prototype is resized to the correct size for \em grid. In case the PDE to be solved is not scalar but a system of equation, the user has to pass the number of equations (\em system_size) to ensure that \em stiffness_matrix_prototype is sized correctly.  
  */
  void ellipse_grid(float_t a, float_t b, std::size_t n_rings, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size = 1);
  
    /** \ingroup fem
      <tt>\#include <fem/utilities.hpp></tt>
      
      Triangulates an arbitrary, planar shape given by \em shape_discretizer and computes a \em grid from the triangulation. The maximum area of each triangle will not exceed \em max_triangle_area.
      In addition, \em stiffness_matrix_prototype is resized to the correct size for \em grid. In case the PDE to be solved is not scalar but a system of equation, the user has to pass the number of equations (\em system_size) to ensure that \em stiffness_matrix_prototype is sized correctly.  
  */
  void triangulate_shape(const BoundaryDiscretizer<2> & shape_discretizer, float_t max_triangle_area, Grid<fem_2d_triangle_types> & grid, ublas::compressed_matrix<float_t> & stiffness_matrix_prototype, std::size_t system_size = 1);
}

#endif

