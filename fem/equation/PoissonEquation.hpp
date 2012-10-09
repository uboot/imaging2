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

#ifndef EQUATION_POISSONEQUATION_H
#define EQUATION_POISSONEQUATION_H

#include <fem/equation/SimpleEquationInterface.hpp>

namespace imaging
{
  /** \ingroup fem_equation
      \brief Implements the Poisson equation with Dirichlet boundary conditions.
      
      This class implements the equation
      \f[
        \nabla\cdot(\nabla u) = f\, , \quad u = v \ \textrm{on}\ \partial \Omega.
      \f]
      Here the force \f$f\f$ is a function on the problem domain
      and the boundary condition \f$v\f$ a function on the boundary of the domain. The user must provide them to assembly the system of linear equations which can be solved for \f$u\f$.
  */
  template<class fem_types>
  class PoissonEquation : public SimpleEquationInterface<fem_types>
  {
    boost::shared_ptr< ublas::vector<float_t> > _rhs_ptr;
    boost::shared_ptr< ublas::mapped_vector<float_t> > _boundary_data_ptr;

  public:
    typedef ublas::fixed_matrix<float_t, fem_types::data_dimension, fem_types::data_dimension> matrix_coefficient_t;

    /** Returns a reference to the vector containing \f$f\f$. */
    const ublas::vector<float_t> & force() const { return *_rhs_ptr; }
    
    /** Sets the vector containing \f$f\f$. */
    void set_force(boost::shared_ptr< ublas::vector<float_t> > rhs_ptr) { _rhs_ptr = rhs_ptr; }
    
    /** Returns a reference to the vector containing \f$v\f$. */
    const ublas::mapped_vector<float_t> & boundary_data() const { return *_boundary_data_ptr; }
    
    /** Sets the vector containing \f$v\f$. */
    void set_boundary_data(boost::shared_ptr< ublas::mapped_vector<float_t> > boundary_data_ptr) { _boundary_data_ptr = boundary_data_ptr; }

    static const size_t boundary_data_type = SimpleEquationInterface<fem_types>::DIRICHLET_DATA;

    static const bool a_active = false;
    static const bool b_active = false;
    static const bool c_active = false;
    static const bool f_active = true;
    static const bool g_active = false;
    
    void stiffness_matrix(std::size_t integrator_node,
                          const FemKernel<fem_types> & kernel,
                          matrix_coefficient_t & A,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & a,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & b,
                          float_t & c) const
    { 
      A = ublas::identity_matrix<float_t>(fem_types::data_dimension);
    }
    
    void force_vector(std::size_t integrator_node,
                      const FemKernel<fem_types> & kernel,
                      float_t & f,
                      ublas::fixed_vector<float_t, fem_types::data_dimension> & g) const
    {
      kernel.grid().interpolate_value(integrator_node, * _rhs_ptr, kernel, f);
    }
    
    void force_vector_at_boundary (std::size_t integrator_node,
                                   const FemKernel< fem_types > & kernel,
                                   float_t & v) const
    {
      kernel.grid().interpolate_boundary_value(integrator_node, * _boundary_data_ptr, kernel, v);
    }
    
    bool sanity_check_stiffness_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      return true;
    }
    
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      if(! (_rhs_ptr && _boundary_data_ptr))
        return true;
        
      if(_rhs_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      if(_boundary_data_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      return true;
    }
    
  };
}


#endif
