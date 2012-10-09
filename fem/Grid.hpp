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

#ifndef GRID_H
#define GRID_H

#include <vector>
#include <map>
#include <set>

#include <fem/FemKernel.hpp>
#include <fem/fem_3d_tetrahedra_types.hpp>
#include <fem/fem_2d_triangle_types.hpp>
#include <fem/fem_1d_types.hpp>
#include <shape/BoundaryDiscretizer.hpp>

// #include <core/cio.hpp>


namespace imaging
{
  template<class fem_types>
  class Image2Grid;
  
  // TODO: clean up accessor functions: void set_vertex() instead of vertex_t & vertex()
  
  /** \ingroup fem
      \brief Provides the data structure and accessor functions for a FE grid.
      
      Grid stores all the data to access the elements of a FE grid and obtain their geometry. In addition it keeps track of the boundary of the grid and elements there.
      
       In the following we give a more detailed description of how Grid stores this information. The basic building block of a grid is an element. We assume an element to have vertices, faces and nodes. Each of these has a unique index (running from zero to Grid::n_element_vertices, from zero to Grid::n_element_faces and from zero to Grid::n_element_nodes) on the reference element which carries over to the transformed elements in the grid.
       
       An element in the grid is identified by its element index which runs from zero to n_elements(). Every element has Grid::n_element_vertices vertices which are identified by their vertex indices. Each of the vertex indices is mapped to a global vertex index which refers to its position in a list of vertex coordinates.In other words, a combination of an <em>element index + vertex index</em> can be mapped to the <em>global vertex index</em> which again identifies the position of the element vertex coordinates in the vertex coordinate list. The vertex coordinates determine the position and the geometry of the element. In a grid consisting of triangular elements an element will always have three vertices. Elements in a quadrilateral or tetrahedral grid are defined by 4 vertices respectively. This also means that the number of vertices per element and the dimension of the vertex coordinates (which are determined by Grid::n_element_vertices and Grid::data_dimension) always have to correspond to the element transformation (Grid::transform_t).
       
       The <em>element nodes</em>, in contrast, are defined by the shape function on the element. Again each node of a element is determined by its <em>node index</em> and is further mapped to a <em>global node index</em>. The coordinates of the nodes do not have to stored because the are completely determined by the coordinates of the vertices of their element. Instead, the global node indices refer to the position of the node in the stiffness matrix and the force vector of the FE problem. For linear shape functions on triangle elements or bilinear shape functions on quadrilateral elements the number of nodes per element is the same as the number of vertices and their actual positions coincide. For higher order shape functions there maybe more nodes than vertices. Clearly, the number of nodes per element (determined by Grid::n_element_nodes) has to match the type of shape functions (Grid::shape_function_t).
       
       Finally, Grid keeps track of the boundary of the FE domain. The term <em>boundary element</em> refers to a face of an element (the <em>parent element</em>) which lies a the boundary of the grid. It is identified by its element index which runs from zero to n_boundary_elements(). For each boundary element the index of its parent element and the index of the corresponding face on the boundary element is stored. A <em>boundary vertex</em> is an element vertex which happens to on a boundary element. The same holds for a <em>boundary node</em>. For each boundary node the outer unit normal at the grid boundary is stored. This information has to be provided during the construction of the grid. 
  */
  template<class fem_types>
  class Grid
  {
  private:
    class BoundaryElement
    {
    public:
      size_t _parent_element;
      size_t _parent_element_face;
      ublas::fixed_vector<float_t, fem_types::data_dimension> _normal;
    };
    
  public:
    /** \brief The number of coordinates of the element vertices (determined by the template parameter \em fem_types).*/
    static const size_t data_dimension = fem_types::data_dimension;
    
    /** \brief The number of vertices per element (determined by the template parameter \em fem_types). */
    static const size_t n_element_vertices = fem_types::transform_t::n_element_vertices;
    
    /** \brief The number of nodes per element (determined by the template parameter \em fem_types). */
    static const size_t n_element_nodes = fem_types::shape_function_t::n_element_nodes;
    
    /** \brief The number of nodes per element (determined by the template parameter \em fem_types). */
    static const size_t n_element_faces = fem_types::transform_t::n_element_faces;
    
    /** \brief The data type of the element coordinates (determined by the template parameter \em fem_types). */
    typedef ublas::fixed_vector<float_t, fem_types::data_dimension> vertex_t;
    
    /** \brief The element transformation class (determined by the template parameter \em fem_types). */
    typedef typename fem_types::transform_t transform_t;
    
    /** \brief The shape function class (determined by the template parameter \em fem_types). */
    typedef typename fem_types::shape_function_t shape_function_t;
    
    /** \brief The integrator class (determined by the template parameter \em fem_types). */
    typedef typename fem_types::integrator_t integrator_t;
    
    /** \brief The boundary integrator class (determined by the template parameter \em fem_types). */
    typedef typename fem_types::boundary_integrator_t boundary_integrator_t;
  
  private:
    typedef ublas::fixed_vector<size_t, n_element_vertices> element_vertices_t;
    typedef ublas::fixed_vector<size_t, n_element_nodes> element_nodes_t;
    
    std::vector<vertex_t> _vertices;
    std::vector<element_vertices_t> _element_vertices;
    std::vector<element_nodes_t> _element_nodes;
    std::vector<BoundaryElement> _boundary_elements;
    std::multimap<size_t, size_t> _boundary_node_boundary_elements;
    std::set<size_t> _boundary_nodes;

    size_t _n_nodes;

    bool _is_regular;

  public:
    
    /** Default constructor. */
    Grid() : _n_nodes(0), _is_regular(false) {}

    /** Returns the number of elements of the grid (excluding boundary elements). */
    size_t n_elements() const { return _element_vertices.size(); }
    
    /** Returns the number of element vertices of the grid. In case there is exactly shape function one per element vertex this equals the number of nodes. */
    size_t n_vertices() const { return _vertices.size(); }
    
    /** Returns the number of element nodes of this grid. The number of element nodes equals the number of degrees of freedom of FE discretization of a scalar PDE on this grid. For a system of equations this number has to be multiplied by the number of equations to obtain the total number of degrees of freedom. */
    size_t n_nodes() const { return _n_nodes; }
    
    /** Returns the number of nodes which are at the boundary of the grid. The equations for these nodes will be influenced by boundary conditions. */
    size_t n_boundary_nodes() const { return _boundary_nodes.size(); }
    
    /** Returns the number of boundary elements, i.e. those boundary elements of parent elements which are at the actual boundary of the grid. */
    size_t n_boundary_elements() const { return _boundary_elements.size(); }
    
    /** Sets the vertex with global index \em global_vertex_index to \em vertex. */
    void set_vertex(size_t global_vertex_index, const vertex_t & vertex) { _vertices[global_vertex_index] = vertex; }
    
    /** Returns a const reference to the vertex with global index \em global_vertex_index. */
    const vertex_t & vertex(size_t global_vertex_index) const { return _vertices[global_vertex_index]; }
    
    /** Sets the vertex with index \em vertex_index on the element \em element_index to \em vertex. */
    void set_vertex(size_t element_index, size_t vertex_index, const vertex_t & vertex) { _vertices[global_vertex_index(element_index, vertex_index)] = vertex; }
    
    /** Returns a const reference to the vertex with index \em vertex_index on the element \em element_index. */
    const vertex_t & vertex(size_t element_index, size_t vertex_index) const { return _vertices[global_vertex_index(element_index, vertex_index)]; }
    
    /** Returns the global vertex index of the node with index \em vertex_index on the element \em element_index. This index corresponds to the index of the vertex coordinates in the vertex coordinate list stored by the Grid object. */
    size_t global_vertex_index(size_t element_index, size_t vertex_index) const { return _element_vertices[element_index](vertex_index); }
    
    /** Sets the global vertex index of the vertex with index \em vertex_index on the element \em element_index. This index corresponds to the index of the vertex coordinates in the vertex coordinate list stored by the Grid object. */
    void set_global_vertex_index(size_t element_index, size_t vertex_index, size_t global_vertex_index) { _element_vertices[element_index](vertex_index) = global_vertex_index; }
    
    /** Returns the global node index of the node with index \em node_index on the element \em element_index. This index corresponds to the position of the node in the stiffness matrix and the force vector of the associated FE problem. */
    size_t global_node_index(size_t element_index, size_t node_index) const { return _element_nodes[element_index](node_index); }
    
    /** Sets the global node index of the node with index \em node_index on the element \em element_index. This index corresponds to the position of the node in the stiffness matrix and the force vector of the associated FE problem. */
    void set_global_node_index(size_t element_index, size_t node_index, size_t global_node_index) { _element_nodes[element_index](node_index) = global_node_index; }
    
    /** Sets the boundary element \em boundary_element_index. Note that the nodes on the boundary element are \em not automatically set to be boundary nodes; the user has to do this manually during grid construction. */
    void set_boundary_element(size_t boundary_element_index, size_t parent_element_index, size_t parent_element_face)
    { 
      transform_t transform;
      ublas::fixed_vector<float_t, fem_types::data_dimension - 1> in;
      vertex_t out;
      element_transform(parent_element_index, transform);
      transform.boundary_normal(parent_element_face, out);
      
      _boundary_elements[boundary_element_index]._parent_element = parent_element_index; 
      _boundary_elements[boundary_element_index]._parent_element_face = parent_element_face;
      _boundary_elements[boundary_element_index]._normal = out;
      
      for(size_t i = 0; i < shape_function_t::n_shape_face_nodes; ++i)
      {
        _boundary_nodes.insert(global_node_index(parent_element_index, shape_function_t::face_node(parent_element_face, i)));
        _boundary_node_boundary_elements.insert(std::pair<size_t, size_t>(global_node_index(parent_element_index, shape_function_t::face_node(parent_element_face, i)), boundary_element_index));
      }
    }
    
    /** Returns the index of the parent element of the boundary element \em boundary_element_index.*/
    size_t parent_element(size_t boundary_element_index) const { return _boundary_elements[boundary_element_index]._parent_element; }
    
    /** Returns the index of the face of the parent element of the boundary element \em boundary_element_index.*/
    size_t parent_element_face(size_t boundary_element_index) const { return _boundary_elements[boundary_element_index]._parent_element_face; }
    
    /** Returns true if the node \em global_node_index lies at the boundary of the grid. */
    bool is_boundary_node(size_t global_node_index) const
    {
      return (_boundary_nodes.count(global_node_index) > 0);
    }
    
    /** Returns true if the node \em node_index on the element \em element_index lies at the boundary of the grid. */
    bool is_boundary_node(size_t element_index, size_t node_index) const
    {
      return is_boundary_node(global_node_index(element_index, node_index));
    }
    
    /** Returns the unit boundary normal of the node \em global_node_index. If the node happens not to be a boundary node an Exception is thrown. */
    const vertex_t boundary_normal(size_t global_node_index) const
    { 
      if( ! _boundary_node_boundary_elements.count(global_node_index))
        throw Exception("Invalid argument 'node' in Grid::boundary_normal().");
        
      std::pair< std::multimap<size_t, size_t>::const_iterator, std::multimap<size_t, size_t>::const_iterator > boundary_elements = _boundary_node_boundary_elements.equal_range(global_node_index);
      
      std::multimap<size_t, size_t>::const_iterator iter = boundary_elements.first;
      
      vertex_t result(0.0);
      for( ; iter != boundary_elements.second; ++iter)
        result += _boundary_elements[iter->second]._normal;
      
      return result / norm_2(result);
    }
  
    /** Returns the unit boundary normal of the node \em node_index on the element \em element_index. If the node happens not to be a boundary node an Exception is thrown. */
    const vertex_t boundary_normal(size_t element_index, size_t node_index) const
    { 
      return boundary_normal(global_node_index(element_index, node_index));
    }

    /** Returns true if the grid is regular, i.e. the geometry of each of its elements is the identical modulo rigid transformation. Marking a grid as regular speeds up the assembly of the stiffness matrix and the force vector. Prominent examples of regular grids are pixel and voxel discretizations. */
    void set_regular(bool is_regular) { _is_regular = is_regular; }
    
    /** Returns true if the grid is regular, i.e. the geometry of each of its elements is the identical modulo rigid transformation. Marking a grid as regular speeds up the assembly of the stiffness matrix and the force vector. Prominent examples of regular grids are pixel and voxel discretizations. */
    bool is_regular() const { return _is_regular; }

    /** Sets the dimensions of the grid. The parameters refer to the \em total numbers of vertices, elements, boundary elements and nodes, respectively. Note that vertices or nodes which belong to more than one element are only counted once. */
    void set_dimensions(size_t n_vertices, size_t n_elements, size_t n_boundary_elements, size_t n_nodes)
    {
      _vertices.resize(n_vertices);
      _element_vertices.resize(n_elements);
      _element_nodes.resize(n_elements);
      _boundary_elements.resize(n_boundary_elements);

      _n_nodes = n_nodes;
    }

    /** Initializes the element transformation \em transform to the element \em element_index. */
    void element_transform(size_t element_index, transform_t & transform) const
    {
      for(size_t i = 0; i < n_element_vertices; ++i)
        transform.assign(i, _vertices[global_vertex_index(element_index, i)]);
    }
    
    /** Interpolates \em data in \em integrator_node on the current element of \em kernel. The result is written to \em value. 
        \param[in] data A vector of size n_nodes(). Its values can be interpreted as the evaluation of a scalar function on the nodes of the grid. interpolate_value() interpolates these values and evaluates the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid element of the grid. */
    void interpolate_value(size_t integrator_node, const ublas::vector<float_t> & data, const FemKernel<fem_types> & kernel, float_t & value) const
    { 
      value = 0.0;
      
      for(size_t i = 0; i < n_element_nodes; ++i)
        value += data(global_node_index(kernel.current_element(), i)) * kernel.shape_value(integrator_node, i);
    }
    
    /** Interpretes \em data as a vector of <em>N</em>-dimensional vectors and interpolates it in \em integrator_node on the current element of \em kernel. The result is written to \em value. 
        \param[in] data A vector of size (<em>N</em> * n_nodes()). Its values can be interpreted as the evaluation of an <em>N</em>-dimensional function on the nodes of the grid. interpolate_value() interpolates these values and evaluates the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid element of the grid. */
    template <size_t N>
    void interpolate_value(size_t integrator_node, const ublas::vector<float_t> & data, const FemKernel<fem_types> & kernel, ublas::fixed_vector<float_t, N> & value) const
    { 
      value.assign(0.0);
      
      for(size_t i = 0; i < n_element_nodes; ++i)
        for(size_t j = 0; j < N; ++j)
          value(j) += data(global_node_index(kernel.current_element(), i) * N + j) * kernel.shape_value(integrator_node, i);
    }
    
    /** Interpolates the gradient of \em data in \em integrator_node on the current element of \em kernel. The result is written to \em value. 
        \param[in] data A vector of size n_nodes(). Its values can be interpreted as the evaluation of a scalar function on the nodes of the grid. interpolate_value() interpolates these values and evaluates the gradient of the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid element of the grid. */
    void interpolate_gradient(size_t integrator_node, const ublas::vector<float_t> & data, const FemKernel<fem_types> & kernel, ublas::fixed_vector<float_t, data_dimension> & gradient) const
    { 
      gradient.assign(0.0);
      
      for(size_t i = 0; i < n_element_nodes; ++i)
        gradient += data(global_node_index(kernel.current_element(), i)) * kernel.shape_gradient(integrator_node, i);
    }
    
    
    /** Interpretes \em data as a vector of <em>N</em>-dimensional vectors and interpolates its gradient in \em integrator_node on the current element of \em kernel. The result is written to \em value. 
        \param[in] data A vector of size (<em>N</em> * n_nodes()). Its values can be interpreted as the evaluation of an <em>N</em>-dimensional function on the nodes of the grid. interpolate_value() interpolates these values and evaluates the gradient (a matrix) of the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid element of the grid. */
    template <size_t N>
    void interpolate_gradient(size_t integrator_node, const ublas::vector<float_t> & data, const FemKernel<fem_types> & kernel, ublas::fixed_matrix<float_t, data_dimension, N> & gradient) const
    { 
      gradient.assign(0.0);
      
      for(size_t i = 0; i < n_element_nodes; ++i)
        for(size_t j = 0; j < N; ++j)
          ublas::matrix_column< ublas::fixed_matrix<float_t, data_dimension, N> >(gradient, j) +=
            data(global_node_index(kernel.current_element(), i) * N + j) * kernel.shape_gradient(integrator_node, i);
    }
    
    /** Interpolates \em data in \em boundary_integrator_node on the current boundary element of \em kernel. The result is written to \em value. 
        \param[in] data A sparse vector of size n_nodes(). It must contain values at positions which correspond to boundary nodes. Its values can be interpreted as the evaluation of a scalar function on the boundary nodes of the grid. interpolate_value() interpolates these values and evaluates the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid boundary element of the grid. */
    void interpolate_boundary_value(size_t boundary_integrator_node, const ublas::mapped_vector<float_t> & data, const FemKernel<fem_types> & kernel, float_t & value) const
    { 
      value = 0.0;
      size_t parent_element = _boundary_elements[kernel.current_boundary_element()]._parent_element;

      for(size_t i = 0; i < n_element_nodes; ++i)
        value += data(global_node_index(parent_element, i)) * kernel.shape_boundary_value(boundary_integrator_node, i);
    }
    
    /** Interpretes \em data as a vector of <em>N</em>-dimensional vectors and interpolates it in \em boundary_integrator_node on the current boundary element of \em kernel. The result is written to \em value. 
        \param[in] data A sparse vector of size (\em N * n_nodes()). It must contain values at positions which correspond to boundary nodes. Its values can be interpreted as the evaluation of an <em>N</em>-dimensional function on the boundary nodes of the grid. interpolate_value() interpolates these values and evaluates the interpolation. 
        \param[in] kernel A FemKernel object. It must be initialized to a valid boundary element of the grid. */
    template <size_t N>
    void interpolate_boundary_value(size_t boundary_integrator_node, const ublas::mapped_vector<float_t> & data, const FemKernel<fem_types> & kernel, ublas::fixed_vector<float_t, N> & value) const
    { 
      value.assign(0.0);
      size_t parent_element = _boundary_elements[kernel.current_boundary_element()]._parent_element;

      for(size_t i = 0; i < n_element_nodes; ++i)
        for(size_t j = 0; j < N; ++j)
          value(j) += data(global_node_index(parent_element, i) * N + j) * kernel.shape_boundary_value(boundary_integrator_node, i);
    }
    
    /** Interpolates \em the tangential derivative of data in \em boundary_integrator_node on the current boundary element of \em kernel. The result is written to \em value. 
        \param[in] data A sparse vector of size n_nodes(). It must contain values at positions which correspond to boundary nodes. Its values can be interpreted as the evaluation of a scalar function on the boundary nodes of the grid. interpolate_value() interpolates these values and evaluates the tangential derivative of the interpolation. Note that it is not possible to compute the normal derivate if \em data is known only on the boundary.
        \param[in] kernel A FemKernel object. It must be initialized to a valid boundary element of the grid. */
    void interpolate_boundary_derivative(size_t boundary_integrator_node, const ublas::mapped_vector<float_t> & data, const FemKernel<fem_types> & kernel, float_t & value) const
    { 
      value = 0.0;
      size_t parent_element = _boundary_elements[kernel.current_boundary_element()]._parent_element;

      for(size_t i = 0; i < n_element_nodes; ++i)
        value += data(global_node_index(parent_element, i)) * kernel.shape_boundary_derivative(boundary_integrator_node, i);
    }
    
    /** Interpretes \em data as a vector of <em>N</em>-dimensional vectors and interpolates it in \em boundary_integrator_node on the current boundary element of \em kernel. The result is written to \em value. 
        \param[in] data A sparse vector of size (\em N * n_nodes()). It must contain values at positions which correspond to boundary nodes. Its values can be interpreted as the evaluation of a <em>N</em>-dimensional function on the boundary nodes of the grid. interpolate_value() interpolates these values and evaluates the tangential derivative of the interpolation. Note that it is not possible to compute the normal derivate if \em data is known only on the boundary.
        \param[in] kernel A FemKernel object. It must be initialized to a valid boundary element of the grid. */
    template <size_t N>
    void interpolate_boundary_derivative(size_t boundary_integrator_node, const ublas::mapped_vector<float_t> & data, const FemKernel<fem_types> & kernel, ublas::fixed_vector<float_t, N> & value) const
    { 
      value.assign(0.0);
      size_t parent_element = _boundary_elements[kernel.current_boundary_element()]._parent_element;
      
      for(size_t i = 0; i < n_element_nodes; ++i)
        for(size_t j = 0; j < N; ++j)
          value(j) += data(global_node_index(parent_element, i) * N + j) * kernel.shape_boundary_derivative(boundary_integrator_node, i);
    }
    
    /** Computes the global coordinates of the \em integrator_node on the current element of \em kernel. The result is written to \em coordinates. */
    void global_coordinates(size_t integrator_node, const FemKernel<fem_types> & kernel, ublas::fixed_vector<float_t, data_dimension> & coordinates) const
    { 
      transform_t transform;
      element_transform(kernel.current_element(), transform);
      
      transform.value( integrator_t().node(integrator_node), coordinates);
    }
    
//     void print_grid() const
//     {
//       std::cout << _vertices << std::endl;
//       std::cout << _element_vertices << std::endl;
//       std::cout << _element_nodes << std::endl;
//     }
  }
  ;
}


#endif
