// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_TRANSFORM_H
#define FEM_TRANSFORM_H

#include <fem/ShapeFunction.hpp>

namespace imaging
{
  /** \ingroup fem 
      \brief Abstract base class of all transformations of the reference element to an element of the FE grid.
  
      All implementations of element transformations should be derived from this class. The template parameters \em N_VERTICES and \em N refer to the number of vertices of the elements and to the dimension of the element respectively.
      
      The derived classes provide functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      To make this more clear, assume the reference element to be the square [-1, 1] x [-1, 1]. I.e. the grid consists of quadrilateral elements. The boundary reference elements of the reference square are the intervals [-1, 1]. First the corresponding Square2dTransform object must be initialized to the 4 vertices of the element it should transform to by calling assign(). It then transforms coordinates in the reference square to coordinates in the element in the 2-dimensional plane. Call value() to perform this transformation and derivative() for its derivative.. Furthermore, the class transforms coordinates in [-1, 1] (i.e. scalar values) to the reference element. This is done by calling boundary2element(). 
      
      \sa Square2dTransform, Triangle2dTransform, Interval1dTransform
  */
  template<size_t N_VERTICES, size_t N_FACES, size_t N>
  class Transform
  {
  protected:
    /** \brief The vertices of the element. */
    ublas::fixed_vector< ublas::fixed_vector<float_t, N>, N_VERTICES> _vertices;

  public:
    /** The number of vertices of the reference element of this transformation. */
    static const size_t n_element_vertices = N_VERTICES;
    
    /** The number of faces of the reference element of this transformation. */
    static const size_t n_element_faces = N_FACES;
    
    /** Returns the vertex index (on the element) of the vertex determined by \em face_index and \em face_vertex_index (on the boundary reference element). E.g. the vertex 1 (there are only vertices 0 and 1 on the boundary reference element) on the face 2 of the square reference element (there are 4 faces with indices from zero to 3) is mapped to vertex 3 on the reference element. */
    static size_t face_vertex(size_t face_index, size_t face_vertex_index);

    /** Sets the vertex \em vertex_index to \em vertex_coordinates. */
    void assign(size_t vertex_index, const ublas::fixed_vector<float_t, N> & vertex_coordinates)
    { _vertices(vertex_index) = vertex_coordinates; }

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, N> & value(const ublas::fixed_vector<float_t, N> & in, ublas::fixed_vector<float_t, N> & out) const;
    
    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, N, N> & derivative(const ublas::fixed_vector<float_t, N> & in, ublas::fixed_matrix<float_t, N, N> & out) const;
    
    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, N, N - 1> & boundary_derivative(size_t face_index, const ublas::fixed_vector<float_t, N - 1> & in, ublas::fixed_matrix<float_t, N, N - 1> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, N> & boundary_normal(size_t face_index, ublas::fixed_vector<float_t, N> & out) const;
    
    /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, N> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, N - 1> & in, ublas::fixed_vector<float_t, N> & out) const;
  };

  /** \ingroup fem 
      \brief Transformation of the square reference element.  
  
      This class implements the transformation of the square [-1, 1] x [-1, 1] (the <em>reference element</em>) to quadrilateral elements and the transformation of the interval [-1, 1] (the <em>boundary reference element</em>) to the faces of this square.
      
      The class provides functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      \sa Transform, Triangle2dTransform, Interval1dTransform
  */
  class Square2dTransform : public Transform<4, 4, 2>
  {
  
  public: 
    static size_t face_vertex(size_t face, size_t vertex);

    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 2, 1> & boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 2, 1> & out) const;

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, 2> & value(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const;

    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 2, 2> & derivative(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 2, 2> & out) const;
    
    /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, 2> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 2> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, 2> & boundary_normal(size_t face, ublas::fixed_vector<float_t, 2> & out) const;
  };

  /** \ingroup fem 
      \brief Transformation of the triangle reference element.  
  
      This class implements the transformation of the triangle determined by the vertices (0, 0), (1, 0), (0, 1) (the <em>reference element</em>) to triangular elements and the transformation of the interval [-1, 1] (the <em>boundary reference element</em>) to the faces of this triangle.
      
      The class provides functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      \sa Transform, Square2dTransform, Interval1dTransform
  */
  class Triangle2dTransform : public Transform<3, 3, 2>
  {
    public:
    static size_t face_vertex(size_t face, size_t vertex);

    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 2, 1> & boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 2, 1> & out) const;

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, 2> & value(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const;

    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 2, 2> & derivative(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 2, 2> & out) const;
    
    /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, 2> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 2> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, 2> & boundary_normal(size_t face, ublas::fixed_vector<float_t, 2> & out) const;
  };

  /** \ingroup fem 
      \brief Transformation of the tetrahedra reference element.  
  
      This class implements the transformation of the tetrahedra determined by the vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1) (the <em>reference element</em>) to tetrahera elements and the transformation of the interval [-1, 1] (the <em>boundary reference element</em>) to the faces of this tetrahedra.
      
      The class provides functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      \sa Transform, Square2dTransform, Interval1dTransform
  */
  class Tetrahedra3dTransform : public Transform<4, 4, 3>
  {
    static const size_t face_vertex_matrix[4][3];
    
    public:
    static size_t face_vertex(size_t face, size_t vertex);

    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 3, 2> & boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 3, 2> & out) const;

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, 3> & value(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const;

    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 3, 3> & derivative(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_matrix<float_t, 3, 3> & out) const;
    
    /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, 3> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 3> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, 3> & boundary_normal(size_t face, ublas::fixed_vector<float_t, 3> & out) const;
  };


  /** \ingroup fem 
      \brief Transformation of the symmetric unit interval reference element.  
  
      This class implements the transformation of the cube [-1, 1]x[-1, 1]x[-1, 1] (the <em>reference element</em>) to quadratic elements and the transformation of the square (the <em>boundary reference element</em>) to the points of the corners.
      
      The class provides functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. 
	  Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      \sa Transform, Square2dTransform, Triangle2dTransform
  */
   class Cube3dTransform : public Transform<8, 6, 3>
  {
    static const size_t face_vertex_matrix[6][4];
    
    public:
    static size_t face_vertex(size_t face, size_t vertex);

    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 3, 2> & boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 3, 2> & out) const;

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, 3> & value(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const;

    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 3, 3> & derivative(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_matrix<float_t, 3, 3> & out) const;
    
    /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, 3> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 3> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, 3> & boundary_normal(size_t face, ublas::fixed_vector<float_t, 3> & out) const;
  };


  /** \ingroup fem 
      \brief Transformation of the symmetric unit cube reference element. TODO:  
  
      This class implements the transformation of the interval [-1, 1] (the <em>reference element</em>) to interval elements and the transformation of the point set {0} (the <em>boundary reference element</em>) to the end points of these intervals.
      
      The class provides functions to transform coordinates on the reference element to grid coordinates and computes the derivate of this transformation. Furthermore, the classes transform coordinates on the <em>boundary reference elements</em> of the reference element (i.e. the reference elements of the faces of the reference element) to coordinates on the reference element.
      
      \sa Transform, Square2dTransform, Triangle2dTransform
  */
  class Interval1dTransform : public Transform<2, 2, 1>
  {
    public:
    static size_t face_vertex(size_t face, size_t vertex);

    /** Computes the derivate of the element transform along the face \em face_index at \em in (on the boundary reference element) and stores it in \em out. For 1-dimensional elements this derivative is not defined. The function must still defined for technical reasons but will in fact never be called by the FE assembly routines. */
    ublas::fixed_matrix<float_t, 1, 0> & boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 0> & in, ublas::fixed_matrix<float_t, 1, 0> & out) const;

    /** Computes the coordinates of \em in (on the reference element) and stores them in \em out. */
    ublas::fixed_vector<float_t, 1> & value(const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 1> & out) const;

    /** Computes the derivate of the element transform at \em in (on the reference element) and stores it in \em out. */
    ublas::fixed_matrix<float_t, 1, 1> & derivative(const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 1, 1> & out) const;
    
        /** Transforms the coordinates \em in on the boundary reference element of the face \em face_index to coordinates on the reference element and stores them in \em out. */
    ublas::fixed_vector<float_t, 1> & boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 0> & in, ublas::fixed_vector<float_t, 1> & out) const;
    
    /** Computes the unit boundary normal at the face \em face_index at \em in (on the boundary reference element) of the reference element and stores it in \em out. If the boundaries of the element are linear segments (which is probably always the case), then this function does actually not depend on \em in. This means the boundary_normal() merely ressembles a look-up table which matches \em face_index with face normal on the reference element. */
    ublas::fixed_vector<float_t, 1> & boundary_normal(size_t face, ublas::fixed_vector<float_t, 1> & out) const;
  };

}

#endif
