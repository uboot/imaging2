// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEFUNCTIONINTERFACE_H
#define FEM_FEFUNCTIONINTERFACE_H

#include <fem/Grid.hpp>
#include <fem/FemKernel.hpp>


namespace imaging
{
  /** \ingroup fem 
      \brief Abstract base class of FE approximations of functions.
      
      This class defines an interface which must be implemented by classes which provide an FE approximation of a function. By FE approximation, we mean a function which can compute its value for a given integration node on the reference element. In addition to the index of the integration node (which itself is of course not sufficient for a meaningful function evaluation) a FemKernel object, which is initialized to the current element, is passed to the function. You can you use the functions provided by FemKernel to interpolate values known on the element nodes to the values and derivatives in the integration nodes.
  */
  class FeFunctionInterface
  {
    public:
    
    /** The type of the output of the function. This might be a scalar or any other kind for which vector operations are implemented (i.e. the operators +, - and *). */
    typedef float_t value_t;
    
    /** The member \em value must evaluate the function in \em integrator_node. The second argument \em fem_kernel will be initialized to the current element. Thus, all the interpolation functions provided by FemKernel can (and most probably have to) be used to compute the actual value. */ 
    template<class fem_types>
    float_t value(size_t integrator_node, const FemKernel<fem_types> fem_kernel) const;
    
    /** Checks if the dimension of the data corresponds to the dimension of the grid stored in \em kernel.
        
        This function is called by the routines evaluating the function. 
        It should return \em false if the provided data does not match the dimensions of the grid.
        Otherwise, return \em true;
    */ 
    template<class fem_types>
    bool sanity_check(const FemKernel<fem_types> & kernel, std::string & error_message) const { return true; } 
  };
  
  
  /** \ingroup fem 
      <tt>\#include <fem/FeFunctionInterface.hpp></tt>
      
      Integrates \em function over \em grid. The parameter \em function has to implement FeFunctionInterface.
  */
  template<class fem_types, class function_t>
  typename function_t::value_t integrate(const function_t & function, const Grid<fem_types> & grid)
  {
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::transform_t transform_t;

    typedef typename function_t::value_t value_t;
  
    integrator_t integrator;

    FemKernel<fem_types> kernel(grid);
    
    std::string sanity_check_message = "";
    if( ! function.sanity_check(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in integrate() with message '" + sanity_check_message + "'.");

    if(grid.is_regular() && grid.n_elements() > 0)
      kernel.set_element(0);
    
    value_t value = 0.0;
      
    for(std::size_t element = 0; element < grid.n_elements(); ++element)
    {

      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(std::size_t k = 0; k < integrator_t::n_nodes; ++k)
        value += function.value(k, kernel) * 
                 kernel.transform_determinant(k) * 
                 integrator.weight(k);
    }
    
    return value;
  }


 /** \ingroup fem 
      <tt>\#include <fem/FeFunctionInterface.hpp></tt>
      
      Computes the  maximum of \em function over the integrator nodes of \em grid. The parameter \em function has to implement FeFunctionInterface.
  */
  template<class fem_types, class function_t>
  typename function_t::value_t maximum(const function_t & function, const Grid<fem_types> & grid)
  {
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::transform_t transform_t;

    typedef typename function_t::value_t value_t;
  
    integrator_t integrator;
    value_t value = 0.0;

    FemKernel<fem_types> kernel(grid);
    
    std::string sanity_check_message = "";
    if( ! function.sanity_check(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in maximum(() with message '" + sanity_check_message + "'.");

    if(grid.n_elements() > 0)
      kernel.set_element(0);
    
    if(integrator_t::n_nodes > 0)
      value = function.value(0, kernel);
      
    for(std::size_t element = 0; element < grid.n_elements(); ++element)
    {
      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(std::size_t k = 0; k < integrator_t::n_nodes; ++k)
        value = max(value, function.value(k,kernel));
    }
    
    return value;
  }

  /** \ingroup fem 
      <tt>\#include <fem/FeFunctionInterface.hpp></tt>
      
      Computes the  minimum of \em function over the integrator nodes of \em grid. The parameter \em function has to implement FeFunctionInterface.
  */
  template<class fem_types, class function_t>
  typename function_t::value_t minimum(const function_t & function, const Grid<fem_types> & grid)
  {
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::transform_t transform_t;

    typedef typename function_t::value_t value_t;
  
    integrator_t integrator;
    value_t value = 0.0;

    FemKernel<fem_types> kernel(grid);
    
    std::string sanity_check_message = "";
    if( ! function.sanity_check(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in minimum(() with message '" + sanity_check_message + "'.");

    if(grid.n_elements() > 0)
      kernel.set_element(0);
    
    if(integrator_t::n_nodes > 0)
      value = function.value(0, kernel);
      
    for(std::size_t element = 0; element < grid.n_elements(); ++element)
    {
      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(std::size_t k = 0; k < integrator_t::n_nodes; ++k)
        value = min(value, function.value(k,kernel));
    }
    
    return value;
  }



 /** \ingroup fem 
      <tt>\#include <fem/FeFunctionInterface.hpp></tt>
      
      Computes the average value of \em function over \em grid. The parameter \em function has to implement FeFunctionInterface.
  */
  template<class fem_types, class function_t>
  typename function_t::value_t average(const function_t & function, const Grid<fem_types> & grid)
  {   
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::transform_t transform_t;

    typedef typename function_t::value_t value_t;
  
    integrator_t integrator;

    FemKernel<fem_types> kernel(grid);
    
    std::string sanity_check_message = "";
    if( ! function.sanity_check(kernel, sanity_check_message) )
      throw Exception("Exception: sanity check failed in average(() with message '" + sanity_check_message + "'.");

    if(grid.is_regular() && grid.n_elements() > 0)
      kernel.set_element(0);
    
    value_t value = 0.0;
    float_t area = 0.0;
      
    for(std::size_t element = 0; element < grid.n_elements(); ++element)
    {
      if(grid.is_regular())
        kernel.lazy_set_element(element);
      else
        kernel.set_element(element);

      for(std::size_t k = 0; k < integrator_t::n_nodes; ++k)
      {
        value += function.value(k, kernel) * 
                 kernel.transform_determinant(k) * 
                 integrator.weight(k);
        area += kernel.transform_determinant(k) * 
                integrator.weight(k);
      }
    }
    
    return value / area;
  }
}


#endif
