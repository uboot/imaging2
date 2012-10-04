// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEMKERNEL_H
#define FEM_FEMKERNEL_H

#include <core/imaging2.hpp>
#include <core/matrix_utilities.hpp>

namespace imaging
{

  template<class fem_types> class Grid;

  template<std::size_t M, std::size_t N>
  float_t transform_det(const ublas::fixed_matrix<float_t, M, N> & derivative)
  {
    return fabs(determinant(derivative));
  }

  /** \ingroup fem
      \brief Computes the values of the shape functions and the element transformation in the integration nodes. 
      
      FemKernel is a class which is used to evaluate functions in the integration nodes on a given element. The integrator, the shape function and the element type is specified by the template parameter \em fem_types. If a FemKernel is set to an element by calling set_element() or set_boundary_element(), it precomputes the values and derivatives of the shape functions and the element transformation in the integration points. Then you can use the accessor functions to query these values. 
      
      A FemKernel object is automatically created by Assembler::assemble(), Assembler::assemble_stiffness_matrix() or Assembler::assemble_force_vector() and then iteratively set to the elements in a given grid during the assembly process. The initialized kernel is then passed to the equation object of the current FE computation. In the implementation of the equation (cf. EquationInterface) you can query the kernel to compute the equation coefficients  which you then pass back to the assembly function. 
      
    \sa EquationInterface
    */
  template<class fem_types>
  class FemKernel
  {
  private:
    static const std::size_t data_dimension = fem_types::data_dimension;

    static const std::size_t n_element_nodes = fem_types::shape_function_t::n_element_nodes;
    static const std::size_t n_element_vertices = fem_types::transform_t::n_element_vertices;

    typedef typename fem_types::shape_function_t shape_function_t;

    typedef typename fem_types::transform_t transform_t;
    typedef typename fem_types::integrator_t integrator_t;
    typedef typename fem_types::boundary_integrator_t boundary_integrator_t;
    
    ublas::fixed_vector< ublas::fixed_vector<float_t, n_element_nodes>, integrator_t::n_nodes> _shape_values;
    ublas::fixed_vector< ublas::fixed_vector< ublas::fixed_vector<float_t, data_dimension>, n_element_nodes >, integrator_t::n_nodes > _shape_gradients;
    ublas::fixed_vector<float_t, integrator_t::n_nodes> _transform_determinants;

    ublas::fixed_vector< ublas::fixed_vector<float_t, n_element_nodes>, boundary_integrator_t::n_nodes > _shape_boundary_values;
    ublas::fixed_vector< ublas::fixed_vector<float_t, n_element_nodes>, boundary_integrator_t::n_nodes > _shape_boundary_derivatives;
    ublas::fixed_vector<float_t, boundary_integrator_t::n_nodes> _boundary_transform_determinants;
    
    ublas::fixed_vector<float_t, data_dimension> _boundary_normal;
    ublas::fixed_vector< ublas::fixed_vector<float_t, data_dimension>, shape_function_t::n_shape_face_nodes > _boundary_normals_at_nodes;
    

    const Grid<fem_types> & _grid;
    std::size_t _current_element;
    std::size_t _current_boundary_element;

  public:

    /** Construct kernel from \em grid. The kernel object stores only a reference to the grid. When the kernel is later on set to an element by calling void set_element(std::size_t element) or set_boundary_element(std::size_t element) the argument \em element refers to \em grid. */
    FemKernel(const Grid<fem_types> & grid) : _grid(grid) {}

    /** Returns a reference to the grid, which was passed to kernel during construction. The kernel object stores only a reference to the grid. */
    const Grid<fem_types> & grid() const { return _grid; };
    
    /** Initializes the kernel to \em element, computes
          - the values and derivatives of the shape functions and
          - the element transform determinant (Jacobian of the transform),
        and stores them. This function is relatively expensive whereas the retrieval of the computed values later on is cheap. As long as set_element() or lazy_set_element() are not called with a different parameter \em element, current_element() will return \em element.
    */
    void set_element(std::size_t element);
    
    /** Initializes the kernel to \em boundary_element, computes
          - the values and derivatives of the boundary shape functions and
          - the boundary element transform determinant (determinant of the Jacobian of the transform),
        and stores them. This function is relatively expensive whereas the retrieval of the computed values later on is cheap. As long as set_boundary_element() is not called with a different parameter \em element, current_boundary_element() will return \em boundary_element.
    */
    void set_boundary_element(std::size_t boundary_element);

    /** Initializes the kernel to \em element, but does \em not compute new values of the shape and transformation functions. Use this function if you know for sure that the geometry of the new element is exactly the same as the one of the previous element. As long as set_element() or lazy_set_element() are not called with a different parameter \em element, current_element() will return \em element.
    */
    void lazy_set_element(std::size_t element) { _current_element = element; }

    /** Returns the value of the shape function assciated with \em element_node in the integration node \em integrator_node. */
    float_t shape_value(std::size_t integrator_node, std::size_t element_node) const
      { return _shape_values(integrator_node)(element_node); }
    
    /** Returns the gradient of the shape function assciated with \em element_node in the integration node \em integrator_node. */
    const ublas::fixed_vector<float_t, fem_types::data_dimension> & shape_gradient(std::size_t integrator_node, std::size_t element_node) const
      { return _shape_gradients(integrator_node)(element_node); }
    
    /** Returns the value of the transformation determinant (determinant of the Jacobian of the transform) in the integration node \em integrator_node. This value reflects the distortion of the actual element with respect to its reference element. Note that \em integrator_node refers to the corresponding node of the \em boundary integrator. */
    float_t transform_determinant(std::size_t integrator_node) const
      { return _transform_determinants(integrator_node); }
  
    /** Returns the value of shape function associated with \em element_node on the parent element of the current boundary element in the integration node \em integrator_node on the boundary element. This means that the shape function <em>on the element</em> is evaluated <em>at the boundary</em> of the element. Keep in mind that here the element is \em not the current element but the parent element of the current boundary element. I.e. this function depends only on the current boundary element. Note that \em integrator_node refers to the corresponding node of the \em boundary integrator. */  
    float_t shape_boundary_value(std::size_t integrator_node, std::size_t element_node) const
      { return _shape_boundary_values(integrator_node)(element_node); }

    /** Returns the gradient of shape function associated with \em element_node on the parent element of the current boundary element in the integration node \em integrator_node on the boundary element. This means that the gradient of the shape function <em>on the element</em> is evaluated <em>at the boundary</em> of the element. Keep in mind that here the element is \em not the current element but the parent element of the current boundary element. I.e. this function depends only on the current boundary element. Note that \em integrator_node refers to the corresponding node of the \em boundary integrator. */  
    float_t shape_boundary_derivative(std::size_t integrator_node, std::size_t element_node) const
      { return _shape_boundary_derivatives(integrator_node)(element_node); }

    /** Returns the value of the transformation determinant (determinant of the Jacobian of the transform) of the transformation which transforms the current boundary element to the reference element. This value reflects the distortion of the actual boundary element with respect to its reference element. */
    float_t boundary_transform_determinant(std::size_t integrator_node) const
      { return _boundary_transform_determinants(integrator_node); }
      
    /** Returns the boundary normal in the integration node \em integrator_node of the current boundary element. Note that \em integrator_node refers to the corresponding node of the \em boundary integrator. */
    const ublas::fixed_vector<float_t, fem_types::data_dimension> & boundary_normal() const
      { return _boundary_normal; }

    /** Returns the index of the current element. In most situations you should assume that the kernel is initialized to this element. */
    std::size_t current_element() const { return _current_element; }
    
    /** Returns the index of the current boundary element. In most situations you should assume that the kernel is initialized to this boundary element. */
    std::size_t current_boundary_element() const { return _current_boundary_element; }
    
    /** Returns \em true if \em element_node refers to a node index on the reference element which lies at the boundary of the grid in the current element. */
    bool is_boundary_node(std::size_t element_node) const { return _grid.is_boundary_node(_current_element,element_node); }
    
    /** Returns the boundary normal in \em element_node of the current element. The parameter \em element_node refers to the index of the node as in the reference element. */
    const ublas::fixed_vector<float_t, fem_types::data_dimension> boundary_normal_at_node(std::size_t element_node) const { return _grid.boundary_normal(_current_element, element_node); }
  }
  ;

  template <class fem_types>
  void FemKernel<fem_types>::set_element(std::size_t element)
  {
    transform_t transform;
    shape_function_t shape_function;
    integrator_t integrator;

    _current_element = element;

    _grid.element_transform(element, transform);

    for(std::size_t integrator_node = 0; integrator_node < integrator_t::n_nodes; ++integrator_node)
    {
      ublas::fixed_matrix<float_t, data_dimension, data_dimension> derivative, inverse_derivative;
      ublas::fixed_vector<float_t, data_dimension> temp;

      transform.derivative(integrator.node(integrator_node), derivative);
      inverse_derivative = inverse(derivative);

      _transform_determinants(integrator_node) = transform_det(derivative);

      for(std::size_t element_node = 0; element_node < n_element_nodes; ++element_node)
      {
        _shape_values(integrator_node)(element_node) = shape_function.value(element_node, integrator.node(integrator_node));
        shape_function.gradient(element_node, integrator.node(integrator_node), temp);
        _shape_gradients(integrator_node)(element_node) = prod(temp, inverse_derivative);
      }
    }

  }



  template <class fem_types>
  void FemKernel<fem_types>::set_boundary_element(std::size_t element)
  {
    transform_t transform;
    shape_function_t shape_function;
    boundary_integrator_t integrator;
    std::size_t parent_element = _grid.parent_element(element);
    std::size_t parent_element_face = _grid.parent_element_face(element);

    _current_boundary_element = element;

    _grid.element_transform(parent_element, transform);

    transform.boundary_normal(parent_element_face, _boundary_normal);
    
    for(std::size_t integrator_node = 0; integrator_node < boundary_integrator_t::n_nodes; ++integrator_node)
    {
      ublas::fixed_matrix<float_t, data_dimension, data_dimension> derivative, inverse_derivative;
      ublas::fixed_matrix<float_t, data_dimension, data_dimension - 1> boundary_transform_derivative;
      ublas::fixed_vector<float_t, data_dimension> element_coordinates, gradient, boundary_normal, transformed_boundary_normal, temp;
      
      transform.boundary2element(parent_element_face, integrator.node(integrator_node), element_coordinates);

      transform.derivative(element_coordinates, derivative);
      inverse_derivative = inverse(derivative);
      transform.boundary_derivative(parent_element_face, integrator.node(integrator_node), boundary_transform_derivative);
      
      _boundary_transform_determinants(integrator_node) = transform_det(boundary_transform_derivative);

      for(std::size_t element_node = 0; element_node < n_element_nodes; ++element_node)
      {
        _shape_boundary_values(integrator_node)(element_node) =
          shape_function.value(element_node, element_coordinates);

        shape_function.gradient(element_node, element_coordinates, temp);
        gradient = prod(temp, inverse_derivative);

        _shape_boundary_derivatives(integrator_node)(element_node) = inner_prod(gradient, transformed_boundary_normal);
      }
    }
  }
  
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 1, 0> & derivative);

  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 1, 1> & derivative);
   
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 2, 1> & derivative);
   
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 3, 2> & derivative);

}


#endif
