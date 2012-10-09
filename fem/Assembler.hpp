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

#ifndef FEM_ASSEMBLER_H
#define FEM_ASSEMBLER_H

#include <fem/Grid.hpp>
#include <fem/FemKernel.hpp>




namespace imaging
{
  /** \ingroup fem
      \brief Assembles the stiffness matrix and force vector of a FE problem.
      
      The Assembler class provides functions to assemble the stiffness matrix and force vector for a given equation and a given grid.
  */
  class Assembler
  {
    static void clear_matrix(ublas::compressed_matrix<float_t> & matrix);
    
  public:

    /** Assembles the stiffness matrix and the force vector for \em equation on \em grid. This is done in one big loop and thus faster than calling assemble_stiffness_matrix() and assemble_force_vector() separately. For performance reasons the type of \em equation is a template parameter.
    
        The class \em equation_t must implement all the functions defined in EquationInterface. It is also advised to derive \em equation_t from EquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */
    template<class fem_types, class equation_t>
    void assemble(const equation_t & equation, const Grid<fem_types> & grid,
                  ublas::compressed_matrix<float_t> & stiffness_matrix,
                  ublas::vector<float_t> & force_vector) const;
    
    /** Assembles the stiffness matrix for \c equation on \c grid. If you want to assembly both, stiffness matrix and force vector, use assemble() to save computation time. For performance reasons the type of \em equation is a template parameter.
    
        The class \em equation_t must implement all the functions defined in EquationInterface. It is also advised to derive \em equation_t from EquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */            
    template<class fem_types, class equation_t>
    void assemble_stiffness_matrix(const equation_t & equation, const Grid<fem_types> & grid,
                  ublas::compressed_matrix<float_t> & stiffness_matrix) const;
                  
    /** Assembles the force vector for \c equation on \c grid. If you want to assembly both, stiffness matrix and force vector, use assemble() to save computation time. For performance reasons the type of \em equation is a template parameter.
    
        The class \em equation_t must implement all the functions defined in EquationInterface. It is also advised to derive \em equation_t from EquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */
    template<class fem_types, class equation_t>
    void assemble_force_vector(const equation_t & equation, const Grid<fem_types> & grid,
                  ublas::vector<float_t> & force_vector) const;

  }
  ;

  template<class fem_types, class equation_t>
  void Assembler::assemble(const equation_t & equation, const Grid<fem_types> & grid,
                                       ublas::compressed_matrix<float_t> & stiffness_matrix,
                                       ublas::vector<float_t> & force_vector) const
  {
    typedef typename fem_types::shape_function_t shape_function_t;
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::boundary_integrator_t boundary_integrator_t;
    typedef typename fem_types::transform_t transform_t;
    
    if(stiffness_matrix.size1() != equation.system_size() * grid.n_nodes() &&
       stiffness_matrix.size2() != equation.system_size() * grid.n_nodes() )
      throw Exception("Exception: Dimension of stiffness matrix does not agree with grid size in Assembler::assemble()");
      
    clear_matrix(stiffness_matrix);

    force_vector.resize(equation.system_size() * grid.n_nodes(), false);
    force_vector.clear();

    integrator_t integrator;
    boundary_integrator_t boundary_integrator;

    FemKernel<fem_types> kernel(grid);

    if(grid.is_regular() && grid.n_elements() > 0)
      kernel.set_element(0);
      
    std::string sanity_check_message = "";
    if( ! equation.sanity_check_stiffness_matrix(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in Assembler::assemble() with message '" + sanity_check_message + "'.");
    if( ! equation.sanity_check_force_vector(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in Assembler::assemble() with message '" + sanity_check_message + "'.");
      
    for(size_t element = 0; element < grid.n_elements(); ++element)
    {
      float_t value;

      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {
          for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
          {
            for(size_t m = 0; m < equation.system_size(); ++m)
            {
              value = 0.0;
  
              for(size_t k = 0; k < integrator_t::n_nodes; ++k)
                value += equation.stiffness_matrix(l, m, i, j, k, kernel) * 
                         kernel.transform_determinant(k) *
                         integrator.weight(k);
  
              stiffness_matrix(equation.system_size() * grid.global_node_index(element, i) + l,
                               equation.system_size() * grid.global_node_index(element, j) + m) += value;
            }
          }

          value = 0.0;

          for(size_t k = 0; k < integrator_t::n_nodes; ++k)
            value += equation.force_vector(l, i, k, kernel) *
                         kernel.transform_determinant(k) *
                         integrator.weight(k);

          force_vector(equation.system_size() * grid.global_node_index(element, i) + l) += value;
        }
      }
    }


    for(size_t element = 0; element < grid.n_boundary_elements(); ++element)
    {
      kernel.set_boundary_element(element);
      size_t parent_element = grid.parent_element(element);
      
      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {  
          for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
          {
            for(size_t m = 0; m < equation.system_size(); ++m)
            {
              float_t value = 0.0;
  
              for(size_t k = 0; k < boundary_integrator_t::n_nodes; ++k)
                value += equation.stiffness_matrix_at_boundary(l, m, i, j, k, kernel) *
                 kernel.boundary_transform_determinant(k) *
                 boundary_integrator.weight(k);
                
              stiffness_matrix(equation.system_size() * grid.global_node_index(parent_element, i) + l,
                               equation.system_size() * grid.global_node_index(parent_element, j) + m) += value;
            }

          }

          float_t value = 0.0;

          for(size_t k = 0; k < boundary_integrator_t::n_nodes; ++k)
            value += equation.force_vector_at_boundary(l, i, k, kernel) *
                 kernel.boundary_transform_determinant(k) *
                 boundary_integrator.weight(k);

          force_vector(equation.system_size() * grid.global_node_index(parent_element, i) + l) += value;
        }
      }
    }

  }

  template<class fem_types, class equation_t>
  void Assembler::assemble_stiffness_matrix(const equation_t & equation,
                                      const Grid<fem_types> & grid,
                                      ublas::compressed_matrix<float_t> & stiffness_matrix) const
  {
    typedef typename fem_types::shape_function_t shape_function_t;
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::boundary_integrator_t boundary_integrator_t;
    typedef typename fem_types::transform_t transform_t;
   
    if(stiffness_matrix.size1() != equation.system_size() * grid.n_nodes() &&
       stiffness_matrix.size2() != equation.system_size() * grid.n_nodes() )
      throw Exception("Exception: Dimension of stiffness matrix does not agree with grid size in Assembler::assemble()");

    clear_matrix(stiffness_matrix);
    
    integrator_t integrator;
    boundary_integrator_t boundary_integrator;

    FemKernel<fem_types> kernel(grid);

    if(grid.is_regular() && grid.n_elements() > 0)
      kernel.set_element(0);
      
    std::string sanity_check_message = "";
    if( ! equation.sanity_check_stiffness_matrix(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in Assembler::assemble() with message '" + sanity_check_message + "'.");  
      
    for(size_t element = 0; element < grid.n_elements(); ++element)
    {
      float_t value;

      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {
          for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
          {
            for(size_t m = 0; m < equation.system_size(); ++m)
            {
              value = 0.0;
  
              for(size_t k = 0; k < integrator_t::n_nodes; ++k)
                value += equation.stiffness_matrix(l, m, i, j, k, kernel) * 
                         kernel.transform_determinant(k) *
                         integrator.weight(k);
  
              stiffness_matrix(equation.system_size() * grid.global_node_index(element, i) + l,
                               equation.system_size() * grid.global_node_index(element, j) + m) += value;
            }
          }
        }
      }
    }


    for(size_t element = 0; element < grid.n_boundary_elements(); ++element)
    {
      kernel.set_boundary_element(element);
      size_t parent_element = grid.parent_element(element);
      
      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {  
          for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
          {
            for(size_t m = 0; m < equation.system_size(); ++m)
            {
              float_t value = 0.0;
  
              for(size_t k = 0; k < boundary_integrator_t::n_nodes; ++k)
                value += equation.stiffness_matrix_at_boundary(l, m, i, j, k, kernel) *
                 kernel.boundary_transform_determinant(k) *
                 boundary_integrator.weight(k);
                
              stiffness_matrix(equation.system_size() * grid.global_node_index(parent_element, i) + l,
                               equation.system_size() * grid.global_node_index(parent_element, j) + m) += value;
            }
          }
        }
      }
    }

  }

  template<class fem_types, class equation_t>
  void Assembler::assemble_force_vector(const equation_t & equation,
                                      const Grid<fem_types> & grid,
                                      ublas::vector<float_t> & force_vector) const
  {
    typedef typename fem_types::shape_function_t shape_function_t;
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::boundary_integrator_t boundary_integrator_t;
    typedef typename fem_types::transform_t transform_t;
   
    force_vector.resize(equation.system_size() * grid.n_nodes(), false);
    force_vector.clear();

    integrator_t integrator;
    boundary_integrator_t boundary_integrator;

    FemKernel<fem_types> kernel(grid);

    if(grid.is_regular() && grid.n_elements() > 0)
      kernel.set_element(0);
      
    std::string sanity_check_message = "";
    if( ! equation.sanity_check_force_vector(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in Assembler::assemble() with message '" + sanity_check_message + "'.");
      
    for(size_t element = 0; element < grid.n_elements(); ++element)
    {
      float_t value;

      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {
          value = 0.0;

          for(size_t k = 0; k < integrator_t::n_nodes; ++k)
            value += equation.force_vector(l, i, k, kernel) *
                         kernel.transform_determinant(k) *
                         integrator.weight(k);

          force_vector(equation.system_size() * grid.global_node_index(element, i) + l) += value;
        }
      }
    }


    for(size_t element = 0; element < grid.n_boundary_elements(); ++element)
    {
      kernel.set_boundary_element(element);
      size_t parent_element = grid.parent_element(element);
      
      for(size_t i = 0; i < fem_types::shape_function_t::n_element_nodes; ++i)
      {
        for(size_t l = 0; l < equation.system_size(); ++l)
        {  
          float_t value = 0.0;

          for(size_t k = 0; k < boundary_integrator_t::n_nodes; ++k)
            value += equation.force_vector_at_boundary(l, i, k, kernel) *
                 kernel.boundary_transform_determinant(k) *
                 boundary_integrator.weight(k);

          force_vector(equation.system_size() * grid.global_node_index(parent_element, i) + l) += value;
        }
      }
    }

  }
}


#endif
