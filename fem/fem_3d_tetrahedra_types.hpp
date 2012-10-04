// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEM3DTRIANGLETYPES_H
#define FEM_FEM3DTRIANGLETYPES_H

#include <fem/ShapeFunction.hpp>
#include <fem/Transform.hpp>
#include <fem/ElementIntegrator.hpp>


namespace imaging
{

  /** \ingroup fem
      \brief FE types for a 3-dimensional FE problem with tetrahedral elements.
      
      This class provides FE traits for 3-dimensional FE problems on tetrahedral elements. It has the following characteristics:
        - reference element: tetrahedra in the space with vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
        - integrator: barycenter integrator with 1 node,
        - boundary integrator: Gaussian integrator with 1 nodes.
        
      \sa fem_1d_types, fem_2d_square_types , fem_2d_triangle_types, fem_3d_cube_types
  */ 
  class fem_3d_tetrahedra_types
  {
  public:
    static const size_t data_dimension = 3;

    typedef Tetrahedra3dTransform transform_t;
    typedef Linear3dShapeFunction shape_function_t;
    typedef TetrahedraIntegrator<1> integrator_t;
    typedef TriangleIntegrator<1> boundary_integrator_t;
  }
  ;

}

#endif
