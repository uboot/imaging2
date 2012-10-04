// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_FEM1DTYPES_H
#define FEM_FEM1DTYPES_H

#include <fem/ShapeFunction.hpp>
#include <fem/Transform.hpp>
#include <fem/ElementIntegrator.hpp>


namespace imaging
{

  /** \ingroup fem
      \brief FE types for a 1-dimensional FE problem.
      
      This class provides FE types for 1-dimensional FE problems. It has the following characteristics:
        - reference element: the interval [-1, 1],
        - integrator: Gaussian integrator with 2 nodes.
        
      \sa fem_2d_square_types, fem_2d_triangle_types, fem_3d_tetrahedra_types, fem_3d_cube_types
  */    
  class fem_1d_types
  {
  public:
    static const std::size_t data_dimension = 1;

    typedef Interval1dTransform transform_t;
    typedef Linear1dShapeFunction shape_function_t;
    typedef IntervalIntegrator<2> integrator_t;
    typedef PointIntegrator boundary_integrator_t;
  }
  ;

}

#endif
