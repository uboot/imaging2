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

#ifndef EQUATION_SIMPLEEQUATIONINTERFACE_H
#define EQUATION_SIMPLEEQUATIONINTERFACE_H

#include <fem/FemKernel.hpp>

namespace imaging
{
  /** \ingroup fem_equation
      \brief Abstract base class of scalar, elliptic PDEs in divergence form.
      
      Derive this interface to solve a scalar, elliptic PDE in divergence form, i.e.
      \f[
        - \nabla \cdot (A \nabla u + \vec a u ) + \vec b \cdot \nabla u + c u = f - \nabla \cdot \vec g
        \quad \textrm{on} \quad \Omega\,.
      \f]
      The weak formulation of the above problem reads as
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \phi + \int_\Omega \vec b \cdot \nabla u \phi +
        \int_\Omega c u \phi = \int_\Omega f \phi +  \int_\Omega \vec g \nabla \phi \,,
      \f]
      for a test function \f$\phi\f$ \em with compact support.
      

      Depending on the boundary conditions (i.e. on the value of the constant
      \c boundary_data_type ), this class implements five different weak formulations of
      boundary value problems.
      
      - <tt> boundary_data_type == NO_BOUNDARY_DATA </tt>
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \psi + \int_\Omega \vec b \cdot \nabla u \psi +
        \int_\Omega c u \psi = \int_\Omega f \psi +  \int_\Omega \vec g \nabla \psi \,,
      \f]
      for a test function \f$\psi\f$ (\em without compact support).
      This implies the boundary conditions
      \f[
        \frac{\partial}{\partial \vec n}(A \nabla u + \vec a u) = 0 \quad \textrm{on} \quad \partial \Omega\,.
      \f]
      
      - <tt> boundary_data_type == IMPLICIT_NEUMANN_DATA </tt>
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \psi + \int_\Omega \vec b \cdot \nabla u \psi +
        \int_\Omega c u \psi = \int_\Omega f \psi +  \int_\Omega \vec g \nabla \psi - 
        \int_{\partial \Omega} v \psi \,,
      \f]
      for a test function \f$\psi\f$ (\em without compact support).
      This implies the boundary conditions
      \f[
        \frac{\partial}{\partial \vec n}(A \nabla u + \vec a u) = v \quad \textrm{on} \quad \partial \Omega\,.
      \f]
      
      - <tt> boundary_data_type == NEUMANN_DATA </tt>
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \phi + \int_\Omega \vec b \cdot \nabla u \phi +
        \int_\Omega c u \phi = \int_\Omega f \phi +  \int_\Omega \vec g \nabla \phi \,,
      \f]
      for a test function \f$\phi\f$ (\em with compact support)
      and
      \f[
        \int_{\partial \Omega} \frac{\partial u}{\partial \vec n} \psi = 
        \int_{\partial \Omega} v \psi\,,
      \f]
      for a test function \f$\psi\f$ (\em without compact support).
      This implies the boundary conditions
      \f[
         \frac{\partial u}{\partial \vec n} = v \quad \textrm{on} \quad \partial \Omega\,.
      \f]
      
      - <tt> boundary_data_type == DIRICHLET_DATA </tt>
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \phi + \int_\Omega \vec b \cdot \nabla u \phi +
        \int_\Omega c u \phi = \int_\Omega f \phi +  \int_\Omega \vec g \nabla \phi \,,
      \f]
      for a test function \f$\phi\f$ (\em with compact support)
      and
      \f[
        \int_{\partial \Omega} u \psi = 
        \int_{\partial \Omega} v \psi\,,
      \f]
      for a test function \f$\psi\f$ (\em without compact support).
      This implies the boundary conditions
      \f[
         u = v \quad \textrm{on} \quad \partial \Omega\,.
      \f]
      
      - <tt> boundary_data_type == MIXED_DATA </tt>
      \f[
        \int_\Omega (A \nabla u + \vec a u ) \nabla \phi + \int_\Omega \vec b \cdot \nabla u \phi +
        \int_\Omega c u \phi = \int_\Omega f \phi +  \int_\Omega \vec g \nabla \phi \,,
      \f]
      for a test function \f$\phi\f$ (\em with compact support)
      and
      \f[
        \int_{\partial \Omega}  \frac{\partial u}{\partial \vec n} \psi + \int_{\partial \Omega} h u \psi = 
        \int_{\partial \Omega} v \psi\,,
      \f]
      for a test function \f$\psi\f$ (\em without compact support).
      This implies the boundary conditions
      \f[
         \frac{\partial u}{\partial \vec n} + h u = v \quad \textrm{on} \quad \partial \Omega\,.
      \f]
  */
  
  template<class fem_types_t>
  class SimpleEquationInterface
  {
  public:
    /** The \em fem_types for which this class template was instantiated. */
    typedef fem_types_t fem_types;
  
    /** Defines the type of matrix as passed to A(). If A() returns a diagonal matrix choosing an appropriate type to store this matrix might yield a speed-up (in particular in higher dimensions; in the plane the difference will most probably neglible). */
    typedef ublas::fixed_matrix<float_t, fem_types::data_dimension, fem_types::data_dimension> matrix_coefficient_t;

    /** Different kinds of boundary data. A detailed description is given in the general part of the class documentation. */
    enum boundary_data_types {
      NO_BOUNDARY_DATA /** No boundary data. This implies a homogenous condition depending on the equation. */,
      IMPLICIT_NEUMANN_DATA /** The same as \em NO_BOUNDARY_DATA, but with a possible non-homogenous condition. */,
      NEUMANN_DATA /** Neumann conditions. */,
      DIRICHLET_DATA /** Dirichlet conditions. */,
      MIXED_DATA  /** Robin conditions. */
    };

    /** Must be set to \em true if \f$a\f$ can be non-zero. If it is \em false, \f$a\f$ will not be evaluated in the assembly functions to save computation time. */
    static const bool a_active = true;
    
    /** Must be set to \em true if \f$b\f$ can be non-zero. If it is \em false, \f$b\f$ will not be evaluated in the assembly functions to save computation time. */
    static const bool b_active = true;
    
    /** Must be set to \em true if \f$c\f$ can be non-zero. If it is \em false, \f$c\f$ will not be evaluated in the assembly functions to save computation time. */
    static const bool c_active = true;
    
    /** Must be set to \em true if \f$f\f$ can be non-zero. If it is \em false, \f$f\f$ will not be evaluated in the assembly functions to save computation time. */
    static const bool f_active = true;
    
    /** Must be set to \em true if \f$g\f$ can be non-zero. If it is \em false, \f$g\f$ will not be evaluated in the assembly functions to save computation time. */
    static const bool g_active = true;
    
    /** The type of boundary conditions for this equation as defined in SimpleEquationInterface::boundary_data_types. */
    static const size_t boundary_data_type;  
    
    /** Evaluates \f$A\f$, \f$a\f$, \f$b\f$ and \f$c\f$ in \em integrator_node on the current element of \em kernel.
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    void stiffness_matrix(std::size_t integrator_node,
                          const FemKernel<fem_types> & kernel,
                          matrix_coefficient_t & A,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & a,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & b,
                          float_t & c) const;
      
    
    /** Evaluates \f$f\f$ and \f$g\f$ in \em integrator_node on the current element of \em kernel.
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    void force_vector(std::size_t integrator_node,
                      const FemKernel<fem_types> & kernel,
                      float_t & f,
                      ublas::fixed_vector<float_t, fem_types::data_dimension> & g) const; 
    
    /** Evaluates \f$h\f$ in \em integrator_node on the current element of \em kernel.
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    void stiffness_matrix_at_boundary(std::size_t integrator_node,
                                      const FemKernel<fem_types> & kernel,
                                      float_t & h) const {}
    
    
    /** Evaluates \f$v\f$ in \em integrator_node on the current element of \em kernel.
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    void force_vector_at_boundary(std::size_t integrator_node,
                                  const FemKernel<fem_types> & kernel,
                                  float_t & v) const {}
        
    /** Checks if the dimension of the data corresponds to the dimension of the grid stored in \em kernel for the assembly of the stiffness matrix.
        
        This function is called by the assemble routine right before the assembly of the stiffness matrix. 
        It should return \em false if data which is necessary for the assembly of the stiffness matrix is still missing.
        Otherwise, return \em true;
    */ 
    bool sanity_check_stiffness_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const { return true; }  
        
    /** Checks if the dimension of the data corresponds to the dimension of the grid stored in \em kernel for the assembly of the force vector.
        
        This function is called by the assemble routine right before the assembly of the force vector. 
        It should return \em false if data which is necessary for the assembly of the force vector is still missing.
        Otherwise, return \em true;
    */ 
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const { return true; }      
  };

}

#endif
