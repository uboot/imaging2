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

