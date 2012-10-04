#ifndef EQUATION_GEODESICACTIVECONTOURSTEP_H
#define EQUATION_GEODESICACTIVECONTOURSTEP_H

#include <core/utilities.hpp>
#include <fem/equation/SimpleEquationInterface.hpp>
#include <boost/shared_ptr.hpp>

namespace imaging
{
  /** \ingroup fem_equation
      \brief Implements an implicit time step of the geodesic active contour evolution.
      
      This class implements the equation
      \f[
        \frac{u - u^0}\delta = | \nabla u |_\epsilon  \Big( \nabla \cdot \Big(g_f \frac{\nabla u}{| \nabla u |_\epsilon }\Big) + g \nu \Big)\,.
      \f]
      Here \f$u\f$, \f$u^0\f$ and \f$g\f$ are scalar functions on the problem domain and \f$\delta > 0\f$ and
      \f[
         | \nabla u |_\epsilon = \sqrt{|\nabla u| + \epsilon^2}\,.
      \f]
      The <em>edge detector</em> \f$g_f\f$ is defined by
      \f[
         g_f = \exp(- \eta f^2)\,.
      \f]
      The user must provide the initial value \f$u^0\f$, the edge detector \f$g\f$ and the step size \f$\delta\f$ to assemble the system of linear equations which can be solved for \f$u\f$.
      Moreover the values for 
        - the balloon force \f$\nu\f$, 
        - the edge parameter \f$\eta\f$, and 
        - the regularization \f$\epsilon\f$
        
      must be set.
  */
  template<class fem_types>
  class GeodesicActiveContourStep : public SimpleEquationInterface<fem_types>
  {
    boost::shared_ptr< ublas::vector<float_t> > _edges_ptr;
    boost::shared_ptr< ublas::vector<float_t> > _initial_function_ptr;
    float_t _step_size;
    float_t _epsilon;
    float_t _edge_parameter;
    float_t _balloon_force;
    
    float_t regularized_abs(const ublas::fixed_vector<float_t, fem_types::data_dimension > & gradient) const
    {
      return sqrt(square(_epsilon) + inner_prod(gradient, gradient));
    }

  public:
    typedef ublas::fixed_matrix<float_t, fem_types::data_dimension, fem_types::data_dimension> matrix_coefficient_t;

    /** Returns a reference to the vector containing the edge function \f$f\f$. */
    const ublas::vector<float_t> & edges() const { return *_edges_ptr; }
    
    /** Sets the vector containing the edge function \f$f\f$. */
    void set_edges(boost::shared_ptr< ublas::vector<float_t> > edges_ptr) { _edges_ptr = edges_ptr; }
    
    /** Returns a reference to the vector containing the initial function \f$u^0\f$. */
    const ublas::vector<float_t> & initial_function() const { return *_initial_function_ptr; }
    
    /** Sets the vector containing the initial function \f$u^0\f$. */
    void set_initial_function(boost::shared_ptr< ublas::vector<float_t> > initial_function_ptr) { _initial_function_ptr = initial_function_ptr; }
    
    /** Returns the current time step \f$\delta\f$. */
    float_t step_size() const { return _step_size; }
    
    /** Sets the size of the time step \f$\delta\f$. */
    void set_step_size(float_t step_size) { _step_size = step_size; }
    
    /** Returns the regularization \f$\epsilon\f$. */
    float_t epsilon() const { return _epsilon;}
    
    /** Sets the regularization \f$\epsilon\f$. */
    void set_epsilon(float_t epsilon) { _epsilon = epsilon; }
    
    /** Returns the edge parameter \f$\eta\f$. */
    float_t edge_parameter() const { return _edge_parameter;}
    
    /** Sets the edge parameter \f$\eta\f$. */
    void set_edge_parameter(float_t edge_parameter) { _edge_parameter = edge_parameter; }
    
    /** Returns the balloon force \f$\nu\f$. */
    float_t balloon_force() const { return _balloon_force;}
    
    /** Sets the balloon force \f$\nu\f$. */
    void set_balloon_force(float_t balloon_force) { _balloon_force = balloon_force; }
    
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
      float_t local_edge;
      ublas::fixed_vector<float_t, fem_types::data_dimension > local_gradient;
      
      kernel.grid().interpolate_value(integrator_node, *_edges_ptr, kernel, local_edge);
      kernel.grid().interpolate_gradient(integrator_node, *_initial_function_ptr, kernel, local_gradient);
      
      A = _step_size * exp(- _edge_parameter * square(local_edge)) / regularized_abs(local_gradient) *
                ublas::identity_matrix<float_t>(fem_types::data_dimension);
      
      c = 1.0 / regularized_abs(local_gradient);        
      
    }
    
    void force_vector(std::size_t integrator_node,
                      const FemKernel<fem_types> & kernel,
                      float_t & f,
                      ublas::fixed_vector<float_t, fem_types::data_dimension> & g) const
    {
    
      float_t local_edge, local_value;
      ublas::fixed_vector<float_t, fem_types::data_dimension > local_gradient;
      
      kernel.grid().interpolate_value(integrator_node, *_edges_ptr, kernel, local_edge);
      kernel.grid().interpolate_value(integrator_node, *_initial_function_ptr, kernel, local_value);
      kernel.grid().interpolate_gradient(integrator_node, *_initial_function_ptr, kernel, local_gradient);
      
      f =  local_value / regularized_abs(local_gradient) + _step_size * exp(- _edge_parameter * square(local_edge)) * _balloon_force;
    }
    
    bool sanity_check_stiffness_matrix(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      if(! (_edges_ptr && _initial_function_ptr))
        return false;
        
      if(_edges_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      if(_initial_function_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      return true;
    }
    
    bool sanity_check_force_vector(const FemKernel<fem_types> & kernel, std::string & error_message) const
    { 
      if(! (_edges_ptr && _initial_function_ptr))
        return false;
        
      if(_edges_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      if(_initial_function_ptr->size() != kernel.grid().n_nodes())
        return false;
        
      return true;
    }
  };

}


#endif
