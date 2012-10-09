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

#ifndef EQUATION_SIMPLEEQUATIONADAPTOR_H
#define EQUATION_SIMPLEEQUATIONADAPTOR_H


#include <fem/equation/EquationInterface.hpp>
#include <fem/equation/SimpleEquationInterface.hpp>
#include <fem/FemKernel.hpp>


namespace imaging
{
  /** \ingroup fem_equation
      \brief Converts a SimpleEquationInterface object into an EquationInterface object.
      
      This class implements EquationInterface and but can be instantiated for a class implementing SimpleEquationInterface. I.e. it enables SimpleEquationInterface objects to be used in functions requiring EquationInterface objects. In the following example this is illustrated for Assembler::assembly():
  \code
  img::Grid<img::fem_2d_square_types> grid;
  // construct grid...
  
  SomeSimpleEquation<img::fem_2d_square_types> simple_equation;
  // provide data (i.e. boundary data, right hand side) to the equation...
  
  img::Assembler<img::fem_2d_square_types>
  img::ublas::compressed_matrix<img::float_t> stiffness_matrix;
  img::ublas::vector<img::float_t> force_vector;
  
  img::SimpleEquationAdaptor< SomeSimpleEquation<img::fem_2d_square_types> >
    equation(simple_equation);
    
  assembler.assemble(equation, grid, stiffness_matrix, force_vector);
  \endcode
  */
  
  template <class equation_t>
  class SimpleEquationAdaptor : public EquationInterface<typename equation_t::fem_types>
  {
    const equation_t & _equation;
    typedef typename equation_t::fem_types fem_types;
    
  public:
    const std::size_t system_size() const { return 1; }
    
    /** Constructs a SimpleEquationAdaptor from \em equation. The class \em equation_t must implement (and should be derived from) SimpleEquationInterface. */
    SimpleEquationAdaptor(const equation_t & equation) : _equation(equation) {}
    
    float_t stiffness_matrix(std::size_t equation, std::size_t component,
                             std::size_t i, std::size_t j,
                             std::size_t integrator_node,
                             const FemKernel<fem_types> & kernel) const
    {
      float_t value = 0.0;
      
      typename equation_t::matrix_coefficient_t A;
      ublas::fixed_vector<float_t, fem_types::data_dimension> a, b;
      float_t c;
      
      if( ! kernel.is_boundary_node(i) ||
            _equation.boundary_data_type == SimpleEquationInterface<fem_types>::NO_BOUNDARY_DATA ||
            _equation.boundary_data_type == SimpleEquationInterface<fem_types>::IMPLICIT_NEUMANN_DATA )
      {
        _equation.stiffness_matrix(integrator_node, kernel, A, a, b, c);
        
        value += inner_prod( prod(A, kernel.shape_gradient(integrator_node, i)),
                                kernel.shape_gradient(integrator_node, j) );

        if(_equation.a_active)
          value += inner_prod(a, kernel.shape_gradient(integrator_node, j) ) *
                  kernel.shape_value(integrator_node, i);

        if(_equation.b_active)
          value += inner_prod(b, kernel.shape_gradient(integrator_node, i) ) *
                  kernel.shape_value(integrator_node, j);

        if(_equation.c_active)
          value += c *
                   kernel.shape_value(integrator_node, i) *
                   kernel.shape_value(integrator_node, j);
      }
      
      return value;
    }
    
    
    float_t force_vector(std::size_t equation,
                         std::size_t i,
                         std::size_t integrator_node,
                         const FemKernel<fem_types> & kernel) const
    {
      float_t value = 0.0;
      
      float_t f;
      ublas::fixed_vector<float_t, fem_types::data_dimension> g;
      
      if( ! kernel.is_boundary_node(i) ||
            _equation.boundary_data_type == equation_t::NO_BOUNDARY_DATA ||
            _equation.boundary_data_type == equation_t::IMPLICIT_NEUMANN_DATA )
      {
        _equation.force_vector(integrator_node, kernel, f, g);
        
        if(_equation.f_active)
          value += f *
                   kernel.shape_value(integrator_node, i);

        if(_equation.g_active)
          value += inner_prod(g, kernel.shape_gradient(integrator_node, i) );
      }
                
      return value;
    }
    
    
    float_t stiffness_matrix_at_boundary(std::size_t equation, std::size_t component,
                                         std::size_t i, std::size_t j,
                                         std::size_t integrator_node,
                                         const FemKernel<fem_types> & kernel) const
    {
      float_t value = 0.0;
      
      float_t h;
      
      if(_equation.boundary_data_type != equation_t::NO_BOUNDARY_DATA)
      {
        if(_equation.boundary_data_type == equation_t::MIXED_DATA)
        {
          _equation.stiffness_matrix_at_boundary(integrator_node, kernel, h);
          
          value += kernel.shape_boundary_derivative(integrator_node, i) * kernel.shape_boundary_value(integrator_node, j) *
                   kernel.boundary_transform_determinant(integrator_node);

          value += h *
                   kernel.shape_boundary_value(integrator_node, i) * kernel.shape_boundary_value(integrator_node, j);
        }

        if(_equation.boundary_data_type == equation_t::NEUMANN_DATA)
          value += kernel.shape_boundary_derivative(integrator_node, i) * kernel.shape_boundary_value(integrator_node, j);
    
        if(_equation.boundary_data_type == equation_t::DIRICHLET_DATA)
          value += kernel.shape_boundary_value(integrator_node, i) * kernel.shape_boundary_value(integrator_node, j);
      }
                
      return value;
    }
    
    
    float_t force_vector_at_boundary(std::size_t equation,
                                     std::size_t i,
                                     std::size_t integrator_node,
                                     const FemKernel<fem_types> & kernel) const
    {
      float_t value = 0.0;
      
      float_t v;
      
      if(_equation.boundary_data_type != equation_t::NO_BOUNDARY_DATA)
      {
        _equation.force_vector_at_boundary(integrator_node, kernel, v);
        
        if(_equation.boundary_data_type == equation_t::IMPLICIT_NEUMANN_DATA)
          value += -1.0 *
                   v *
                   kernel.shape_boundary_value(integrator_node, i);
        else
          value += v *
                   kernel.shape_boundary_value(integrator_node, i);
      }      
                
      return value;
    }
    
    bool sanity_check_stiffness_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      return _equation.sanity_check_stiffness_matrix(kernel, error_message);
    }
    
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      return _equation.sanity_check_force_vector(kernel, error_message);
    }
  };
}

#endif

