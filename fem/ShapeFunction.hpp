// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H

#include <core/imaging2.hpp>


namespace imaging
{
  /** \ingroup fem
      \brief Abstract base class of all shape functions classes.
      
      All implementations of shape functions on reference elements should be derived from this class. The template parameters \em N and \em N_NODES refer to the dimension of the reference element and to the number of nodes nodes on the element respectively.
      
      \sa Bilinear2dShapeFunction, Linear2dShapeFunction, Linear1dShapeFunction
  */
  template<size_t N_NODES, size_t N_FACE_NODES, size_t N>
  class ShapeFunction
  {
  public:
    /** The number of nodes per element for this shape function. */
    static const size_t n_element_nodes = N_NODES;
    
    /** The number of nodes per face of the reference element of this shape function. */
    static const size_t n_shape_face_nodes = N_FACE_NODES;
    
    /** Returns the node index (on the element) of the node determined by \em face_index and \em face_node_index (on the boundary reference element). E.g. the vertex 1 (there are only vertices 0 and 1 on the boundary reference element) on the face 2 of the square reference element (there are 4 faces with indices from zero to 3) is mapped to vertex 3 on the reference element. */
    static size_t face_node(size_t face, size_t node);
  
    /** Returns the value of the shape function corresponding to \em node_index at \em in. */
    float_t value(size_t node_index, const ublas::fixed_vector<float_t, N> & in) const;
    
    /** Computes the gradient of the shape function corresponding to \em node_index at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, N> & gradient(size_t node_index, const ublas::fixed_vector<float_t, N> & in, ublas::fixed_vector<float_t, N> & out) const;
  };

  /** \ingroup fem
      \brief Bilinear shape function on the square.
      
      Bilinear shape function on the square [-1, 1] x [-1, 1] (i.e. 4 element nodes) You will rarely use this class directly, but might plug it into your own FE traits class as \em shape_function_t.
  */
  class Bilinear2dShapeFunction : public ShapeFunction<4, 2, 2>
  {
  public:
    static size_t face_node(size_t face, size_t node);
    
    /** Returns the value of the shape function corresponding to \em node_index at \em in. */
    float_t value(size_t node, const ublas::fixed_vector<float_t, 2> & in) const;

    /** Computes the gradient of the shape function corresponding to \em node_index at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, 2> & gradient(size_t node_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const;
  };
  

class Trilinear3dShapeFunction : public ShapeFunction<8, 4, 3>
  {
  public:
    static size_t face_node(size_t face, size_t node);
    
    /** Returns the value of the shape function corresponding to \em node_index at \em in. */
    float_t value(size_t node, const ublas::fixed_vector<float_t, 3> & in) const;

    /** Computes the gradient of the shape function corresponding to \em node_index at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, 3> & gradient(size_t node_index, const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const;
  };


  /** \ingroup fem
      \brief Linear shape function on the triangle.
      
      Linear shape function on the triangle determined by the vertices (0, 0), (1, 0), (0, 1) (i.e. 3 element nodes). You will rarely use this class directly, but might plug it into your own FE traits class as \em shape_function_t.
  */
  class Linear2dShapeFunction : public ShapeFunction<3, 2, 2>
  {
  public:
    static size_t face_node(size_t face, size_t node);
    
    /** Returns the value of the shape function corresponding to \em node_index at \em in. */
    float_t value(size_t node_index, const ublas::fixed_vector<float_t, 2> & in) const;


    /** Computes the gradient of the shape function corresponding to \em node_index at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, 2> & gradient(size_t node_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const;

  };
  
/** \ingroup fem
      \brief Linear shape function on the tetrahedra.
      
      Linear shape function on the tetrahedra determined by the vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1) (i.e. 4 element nodes). You will rarely use this class directly, but might plug it into your own FE traits class as \em shape_function_t.
  */
  class Linear3dShapeFunction : public ShapeFunction<4, 3, 3>
  {
  public:
    static size_t face_node(size_t face, size_t node);
    /** Returns the value of the shape function corresponding to \em node_index at \em in. */
    float_t value(size_t node_index, const ublas::fixed_vector<float_t, 3> & in) const;

    /** Computes the gradient of the shape function corresponding to \em node_index at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, 3> & gradient(size_t node_index, const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const;

  };


  /** \ingroup fem
      \brief Linear shape function on the symmetric unit interval.
      
      Linear shape function on the interval [-1, 1]. (i.e. 2 element nodes). You will rarely use this class directly, but might plug it into your own FE traits class as \em shape_function_t.
  */
  class Linear1dShapeFunction : public ShapeFunction<2, 1, 1>
  {
  public:
    static size_t face_node(size_t face, size_t node);
    
    /** Returns the value of the shape function corresponding to \em node at \em in. */
    float_t value(size_t node, const ublas::fixed_vector<float_t, 1> & in) const;


    /** Computes the gradient of the shape function corresponding to \em node at \em in and stores the result in \em out. */
    ublas::fixed_vector<float_t, 1> & gradient(size_t node, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 1> & out) const;
  };

}


#endif
