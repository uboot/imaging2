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
