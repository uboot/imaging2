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
