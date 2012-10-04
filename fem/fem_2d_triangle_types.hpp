// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEM2DTRIANGLETYPES_H
#define FEM_FEM2DTRIANGLETYPES_H

#include <fem/ShapeFunction.hpp>
#include <fem/Transform.hpp>
#include <fem/ElementIntegrator.hpp>


namespace imaging
{

  /** \ingroup fem
      \brief FE types for a 2-dimensional FE problem with triangular elements.
      
      This class provides FE traits for 2-dimensional FE problems on triangular elements. It has the following characteristics:
        - reference element: triangle in the plane with vertices (0, 0), (1, 0), (0, 1),
        - integrator: barycenter integrator with 1 node,
        - boundary integrator: Gaussian integrator with 2 nodes.
        
      \sa fem_1d_types, fem_2d_square_types, fem_3d_tetrahedra_types, fem_3d_cube_types
  */ 
  class fem_2d_triangle_types
  {
  public:
    static const size_t data_dimension = 2;

    typedef Triangle2dTransform transform_t;
    typedef Linear2dShapeFunction shape_function_t;
    typedef TriangleIntegrator<1> integrator_t;
    typedef IntervalIntegrator<2> boundary_integrator_t;
  }
  ;

}

#endif
