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

#ifndef EQUATION_DIFFUSIONSTEP_H
#define EQUATION_DIFFUSIONSTEP_H

#include <fem/equation/SimpleEquationInterface.hpp>
#include <boost/shared_ptr.hpp>


namespace imaging
{
  /** \ingroup fem_equation
      \brief Implements an implicit time step of an isotropic diffusion equation.
      
      This class implements the equation
      \f[
        \frac{u - u^0}\delta = \nabla \cdot (a \nabla u)\,.
      \f]
      Here \f$u\f$, \f$u^0\f$ and \f$a\f$ are scalar functions on the problem domain and \f$\delta > 0\f$. The user must provide the initial value \f$u^0\f$, the diffusion \f$a\f$ and the step size \f$\delta\f$ to assemble the system of linear equations which can be solved for \f$u\f$.
      */
  template<class fem_types>
  class DiffusionStep : public SimpleEquationInterface<fem_types>
  {
    boost::shared_ptr< ublas::vector<float_t> > _input_ptr;
    boost::shared_ptr< ublas::vector<float_t> > _diffusion_ptr;
    float_t _step_size;

  public:
    typedef ublas::fixed_matrix<float_t, fem_types::data_dimension, fem_types::data_dimension> matrix_coefficient_t;
    typedef typename ublas::fixed_vector<float_t, fem_types::data_dimension> vector_coefficient_t;

    DiffusionStep() : _step_size(0.0) {}

    /** Returns a reference to the vector containing the initial value \f$u^0\f$. */
    const ublas::vector<float_t> & input() const { return *_input_ptr; }
    
    /** Sets the vector containing the initial value \f$u^0\f$. */
    void set_input(boost::shared_ptr< ublas::vector<float_t> > input_ptr) { _input_ptr = input_ptr; }
    
    /** Returns a reference to the vector containing the diffusion \f$a\f$. */
    const ublas::vector<float_t> & diffusion() const { return *_diffusion_ptr; }
    
    /** Sets the vector containing the diffusion \f$a\f$. */
    void set_diffusion(boost::shared_ptr< ublas::vector<float_t> > diffusion_ptr) { _diffusion_ptr = diffusion_ptr; }
    
    /** Returns the current time step \f$\delta\f$. */
    float_t step_size() const { return _step_size; }
    
    /** Sets the size of the time step. */
    void set_step_size(float_t step_size) { _step_size = step_size; }
    
    static const size_t boundary_data_type = SimpleEquationInterface<fem_types>::NO_BOUNDARY_DATA;

    static const bool a_active = false;
    static const bool b_active = false;
    static const bool c_active = true;
    static const bool f_active = true;
    static const bool g_active = false;
    
    void stiffness_matrix(std::size_t integrator_node,
                          const FemKernel<fem_types> & kernel,
                          matrix_coefficient_t & A,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & a,
                          ublas::fixed_vector<float_t, fem_types::data_dimension> & b,
                          float_t & c) const
    { 
      float_t local_diffusion;
      kernel.grid().interpolate_value(integrator_node, *_diffusion_ptr, kernel, local_diffusion);
      A = _step_size * local_diffusion * ublas::identity_matrix<float_t>(fem_types::data_dimension);
      c = 1.0;
    }
    
    void force_vector(std::size_t integrator_node,
                      const FemKernel<fem_types> & kernel,
                      float_t & f,
                      ublas::fixed_vector<float_t, fem_types::data_dimension> & g) const
    {
      kernel.grid().interpolate_value(integrator_node, *_input_ptr, kernel, f);
    }
    
    bool sanity_check_stiffness_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      if(! _diffusion_ptr)
        return false;
        
      if(_diffusion_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      return true;
    }
    
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      if(! _input_ptr)
        return false;
        
      if(_input_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      return true;
    }
  };

}


#endif
