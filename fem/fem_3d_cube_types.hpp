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

#ifndef FEM_FEM3DCUBETYPES_H
#define FEM_FEM3DCUBETYPES_H

#include <fem/ShapeFunction.hpp>
#include <fem/Transform.hpp>
#include <fem/ElementIntegrator.hpp>


namespace imaging
{
  
  /** \ingroup fem
      \brief FE types for a 3-dimensional FE problem with cube elements.
      
      This class provides FE types for 3-dimensional FE problems on cube elements. It has the following characteristics:
        - reference element: the cube [-1, -1, -1] x [1, 1, 1]  (diagonal corners) in the space,
        - integrator: Gaussian integrator with 8 nodes,
        - boundary integrator: Gaussian integrator with 4 nodes.
        
      \sa fem_1d_types, fem_2d_triangle_types, fem_2d_square_types, fem_3d_tetrahedra_types
  */ 
  class fem_3d_cube_types
  {
  public:
    static const size_t data_dimension = 3;

    typedef Cube3dTransform transform_t;
    typedef Trilinear3dShapeFunction shape_function_t;
    typedef CubeIntegrator<8> integrator_t;
    typedef SquareIntegrator<4> boundary_integrator_t;

  }
  ;

}

#endif
