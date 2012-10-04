// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_ELEMENTINTEGRATOR_H
#define FEM_ELEMENTINTEGRATOR_H

#include <core/imaging2.hpp>
#include <core/vector_utilities.hpp>
#include <core/matrix_utilities.hpp>


namespace imaging
{
  /** \ingroup fem
      \brief Abstract base class of all integrator classes.
      
      All implementations of integrators on reference elements should be derived from this class. The template parameters \em N and \em N_NODES refer to the dimension of the reference element and to the number of nodes of the integrator respectively.
      
      \sa SquareIntegrator, TriangleIntegrator, IntervalIntegrator, PointIntegrator, Tetraheda
  */
  template<std::size_t N, std::size_t N_NODES>
  class ElementIntegrator
  {
  protected:
    /** \brief The coordinates of the integration nodes on the reference element. */
    ublas::fixed_vector< ublas::fixed_vector<float_t, N>, N_NODES> _nodes;
    
    /** \brief The corresponding weights. */
    ublas::fixed_vector<float_t, N_NODES> _weights;
  public:
     /** \brief The number of integration nodes. */
    static const std::size_t n_nodes = N_NODES;

    /** Returns the coordinates of the <em>i</em>-th integration node on the reference element. */
    const ublas::fixed_vector<float_t, N> & node(std::size_t i) const { return _nodes(i); }
    
    /** Returns the weight corresponding to the <em>i</em>-th integration. */
    float_t weight(std::size_t i) const { return _weights(i); }
  };
  
  


    /** \ingroup fem
        \brief Gaussian integrator on the square.
      
        Gaussian integrator on the square [-1, 1] x [-1, 1] with \em N nodes. SquareIntegrator<2, 4> is used as element integrator in fem_2d_square_types. SquareIntegrator is currently only implemented for \em N_NODES = 4. You will rarely use this class directly, but might plug it into your own FE traits class as \em integrator_t or \em boundary_integrator_t.
  */
  template<std::size_t N_NODES>
  class SquareIntegrator : public ElementIntegrator<2, N_NODES>
  {
  public:
    SquareIntegrator();
  };

  template<>
  SquareIntegrator<4>::SquareIntegrator();
  

    /** \ingroup fem
        \brief Gaussian integrator on the square.
      
        Gaussian integrator on the cube [-1, 1] x [-1, 1] x [-1,1] with \em N nodes. 
		SquareIntegrator<3, 8> is used as element integrator in fem_3d_cube_types. CubeIntegrator is currently only implemented for \em N_NODES = 4. You will rarely use this class directly, but might plug it into your own FE traits class as \em integrator_t or \em boundary_integrator_t.
  */
  template<std::size_t N_NODES>
  class CubeIntegrator : public ElementIntegrator<3, N_NODES>
  {
  public:
    CubeIntegrator();
  };

  template<>
  CubeIntegrator<8>::CubeIntegrator();
  
  /** \ingroup fem
      \brief Gaussian integrator on the triangle.
      
      Gaussian integrator on the triangle determined by the vertices (0, 0), (1, 0), (0, 1). TriangleIntegrator<1, 1> is used as element integrator in fem_2d_triangle_types. TriangleIntegrator is currently implemented for \em N_NODES = 1, 4, 7. You will rarely use this class directly, but might plug it into your own FE traits class as \em integrator_t or \em boundary_integrator_t.
  */
  template<std::size_t N_NODES>
  class TriangleIntegrator : public ElementIntegrator<2, N_NODES>
  {
  public:
    TriangleIntegrator();
  };

  template<>
  TriangleIntegrator<1>::TriangleIntegrator();

  template<>
  TriangleIntegrator<4>::TriangleIntegrator();
  
  template<>
  TriangleIntegrator<7>::TriangleIntegrator();


 /** \ingroup fem
      \brief Gaussian integrator on the tetrahedron.
      
      Gaussian integrator on the tetrahedron determined by the vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1). TetrahedraIntegrator<1, 1> is used as element integrator in fem_3d_tetrahedron_types. TetrahedraIntegrator is currently implemented for \em N_NODES = 1. You will rarely use this class directly, but might plug it into your own FE traits class as \em integrator_t or \em boundary_integrator_t.
  */
  template<std::size_t N_NODES>
  class TetrahedraIntegrator : public ElementIntegrator<3, N_NODES>
  {
  public:
    TetrahedraIntegrator();
  };

  template<>
  TetrahedraIntegrator<1>::TetrahedraIntegrator();

  
  /** \ingroup fem
      \brief Gaussian integrator on the symmetric unit interval.
      
      Gaussian integrator on the interval [-1, 1]. IntervalIntegrator<1> is used as element integrator in fem_1d_types and as boundary element integrator in fem_2d_square_types and fem_2d_triangle_types. IntervalIntegrator is currently implemented for \em N_NODES = 1, 2. You will rarely use this class directly, but might plug it into your own FE traits class as \em integrator_t or \em boundary_integrator_t.
  */
  template<std::size_t N_NODES>
  class IntervalIntegrator : public ElementIntegrator<1, N_NODES>
  {
  public:
    IntervalIntegrator();
  };

  template<>
  IntervalIntegrator<2>::IntervalIntegrator();

  template<>
  IntervalIntegrator<1>::IntervalIntegrator();
  
  /** \ingroup fem
      \brief (Gaussian) integrator on a point.
      
      This integrator is included for technical reasons and is used as boundary integrator in fem_1d_types. This allows the FE assembly code to generically treat boundary points as boundary elements. You will never use this class directly, but have to plut it into your own FE traits class as \em boundary_integrator_t if you solve a 1-dimensional problem.
  */
  class PointIntegrator : public ElementIntegrator<0, 1>
  {
  public:
    PointIntegrator() : ElementIntegrator<0, 1>() { _weights.assign(1.0); }
  };

}


#endif
