// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEM2DSQUARETYPES_H
#define FEM_FEM2DSQUARETYPES_H

#include <fem/ShapeFunction.hpp>
#include <fem/Transform.hpp>
#include <fem/ElementIntegrator.hpp>


namespace imaging
{
  
  /** \ingroup fem
      \brief FE types for a 2-dimensional FE problem with rectangular elements.
      
      This class provides FE types for 2-dimensional FE problems on rectangular elements. It has the following characteristics:
        - reference element: the square [-1, 1] x [-1, 1] in the plane,
        - integrator: Gaussian integrator with 4 nodes,
        - boundary integrator: Gaussian integrator with 2 nodes.
        
      \sa fem_1d_types, fem_2d_triangle_types, fem_3d_tetrahedra_types, fem_3d_cube_types
  */ 
  class fem_2d_square_types
  {
  public:
    static const size_t data_dimension = 2;

    typedef Square2dTransform transform_t;
    typedef Bilinear2dShapeFunction shape_function_t;
    typedef SquareIntegrator<4> integrator_t;
    typedef IntervalIntegrator<2> boundary_integrator_t;
  }
  ;

}

#endif
