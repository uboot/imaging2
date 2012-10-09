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

#ifndef EQUATION_EQUATIONINTERFACE_H
#define EQUATION_EQUATIONINTERFACE_H

#include <fem/FemKernel.hpp>

namespace imaging
{
  /** \ingroup fem_equation
      \brief Abstract base class of all PDEs to be solved using the finite element module.
      
      To solve a (system of) PDEs with the \em imaging2 finite element module the user must design a class which implements the public members of this interface. It is also recommended to actually derive your class from EquationInterface.
      Then, an object of this class should be passed to Assembler::assemble(), Assembler::assemble_stiffness_matrix() or Assembler::assemble_force_vector() to compute the system of linear equation corresponding to the PDE. 
      
      The interface consists of four functions computing the coefficients in the stiffness matrix and the force vector respectively for a given element (determined by the current state of FemKernel). The function system_size() returns the number of equations in the system of PDEs.
      
      For scalar PDEs it is advisable to make use of the interface SimpleEquationInterface. A simple equation can be transformed into a equation complying with EquationInterface by using SimpleEquationAdaptor.
      
      To understand the parameters of the member function we give a (very technical) definition of the PDE problem this interface can model. Imagine that we are interested in the problem
      \f[
      \vec F(\vec u) = \vec f\,,
      \f]
      where \f$\vec F\f$ depends on the vector-valued function \f$\vec u\f$ and its (higher order) derivative. We will rewrite this equation as
      \f[
      F^k\big((u^\ell)_\ell\big) = f^k\,,
      \f]
      where \f$k\f$ and \f$\ell\f$ are indices in the range of the number of equations of the system. 
      Furthermore we assume a basis \f$(\psi_i)_i\f$ of shape functions on the problem domain, where each \f$\psi_i\f$ corresponds to a node on the grid. We define the vector-valued shape functions \f$\vec \psi_i^k\f$ as
      \f[
      \big(\vec \psi_i^k\big)^m = \delta_{km} \psi_i\,.
      \f]
      
      Following the common finite element approach we we express the unknown function \f$\vec u\f$ in terms of the shape functions, i.e.
      \f[
      u^\ell = \sum_i u_i^\ell \psi_i\,.
      \f] 
      We hit the above equation with each of the test functions and integrate the resulting equations:
      \f[
      \int_\Omega \vec F(\vec u) \cdot \vec \psi_j^k = \int_\Omega \vec f  \cdot \vec \psi_j^k\,.
      \f]
      Then we can write the above equation in the following form: 
      \f[
      \sum_i \sum_\ell u_i^\ell \int_\Omega F^{k\ell}_i (\vec u) \psi_j =
      \int_\Omega f^k \psi_j\,.
      \f]
      Applying the Trace Theorem yields an extended version of this equation (to save symbols we re-use \f$F\f$ and \f$f\f$; they are \em not the same functions as above):
      \f[
      \sum_i \sum_\ell u_i^\ell \Big(\underbrace{\int_\Omega F^{k\ell}_i(\vec u, \psi_j)}_{A^{k\ell}_{ij}} +
      \underbrace{\int_{\partial \Omega} G^{k\ell}_i(\vec u, \psi_j)}_{B^{k\ell}_{ij}} \Big) =
      \underbrace{\int_\Omega f^k(\psi_j)}_{a^k_j} + 
      \underbrace{\int_{\partial \Omega} g^k(\psi_j)}_{b^k_j} \,.
      \f]
      The parts \f$A^{k\ell}_{ij}\f$, \f$B^{k\ell}_{ij}\f$, \f$a^k_j\f$ and \f$b^k_j\f$ then correspond to the four members we mentioned in the introduction. Note that the above version of the equation is not uniquely determined. It depends on how the user applies the Trace Theorem, i.e. which parts of the equation you chose to transport to the boundary.
      
      If she also must incorporate boundary conditions she might choose not to evaluate the equation on the boundary at all. This is equivalent to hit the equation with test functions with compact support only. Then, the terms  \f$B^{k\ell}_{ij}\f$ and \f$b^k_j\f$ vanish at first.
      Instead, consider boundary conditions
       \f[
      \vec G(\vec u) = \vec g\quad \textrm{on }\partial \Omega\,.
      \f]
      and hit them with non-compactly supported test functions of the above form, i.e. test functions \f$\vec \psi_i^k\f$ where \f$i\f$ refers to a boundary node.
      Performing the same procedure as above on this equation adds the terms \f$B^{k\ell}_{ij}\f$ and \f$b^k_j\f$ again (this time as a result of the boundary equation and not of the Trace Theorem).
  */
  template <class fem_types>
  class EquationInterface
  {
    
  public:
    /** Returns the number of equations in this system of PDEs. For scalar PDEs this is 1. */
    std::size_t system_size() const;
    
    /** Computes the contribution of the integral \f$A^{k\ell}_{ij}\f$ in \em integrator_node on the current element of \em kernel. The node indices \em i and \em j are node indices on the element, not global indices. 
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via the indices \em i, \em j and \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    float_t stiffness_matrix(std::size_t k, std::size_t l,
                             std::size_t i, std::size_t j,
                             std::size_t integrator_node,
                             const FemKernel<fem_types> & kernel) const;
    
    /** Computes the contribution of the integral \f$a^k_j\f$ in \em integrator_node on the current element of \em kernel. The node index \em i is the node index on the element, not a global index. 
        
        As \em kernel is set to the current element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via the indices \em i and \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    float_t force_vector(std::size_t k,
                         std::size_t i,
                         std::size_t integrator_node,
                         const FemKernel<fem_types> & kernel) const;
  
    /** Computes the contribution of the integral \f$B^{k\ell}_{ij}\f$ in \em integrator_node on the current boundary element of \em kernel. The node indices \em i and \em j are node indices on the element, not global indices. 
        
        As \em kernel is set to the current boundary element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via the indices \em i, \em j and \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    float_t stiffness_matrix_at_boundary(std::size_t k, std::size_t l,
                                         std::size_t i, std::size_t j,
                                         std::size_t integrator_node,
                                         const FemKernel<fem_types> & kernel) const;
        
    /** Computes the contribution of the integral \f$b^k_j\f$ in \em integrator_node on the current element of \em kernel. The node index \em i is the node index on the element, not a global index. 
        
        As \em kernel is set to the current boundary element, in the implementation of this function all relevant values of the shape functions and the element transform can be retrieved from \em kernel via the indices \em i and \em integrator_node. One can also obtain the FE grid from the kernel (FemKernel::grid()) and use its interpolation methods (passing the already correctly initialized \em kernel to them).
    */
    float_t force_vector_at_boundary(std::size_t k,
                                     std::size_t i,
                                     std::size_t integrator_node,
                                     const FemKernel<fem_types> & kernel) const;
        
    /** Checks if the dimension of the data corresponds to the dimension of the grid stored in \em kernel for the assembly of the stiffness matrix.
        
        This function is called by the assemble routine right before the assembly of the stiffness matrix. 
        It should return \em false if data which is necessary for the assembly of the stiffness matrix is still missing.
        Otherwise, return \em true;
    */ 
    bool sanity_check_stiffnesss_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const { return true; }  
        
    /** Checks if the dimension of the data corresponds to the dimension of the grid stored in \em kernel for the assembly of the force vector.
        
        This function is called by the assemble routine right before the assembly of the force vector. 
        It should return \em false if data which is necessary for the assembly of the force vector is still missing.
        Otherwise, return \em true;
    */ 
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const { return true; }                                                            
  };
}

#endif
