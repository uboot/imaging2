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

#ifndef SHAPE_DISCRETIZABLESHAPEINTERFACE_H
#define SHAPE_DISCRETIZABLESHAPEINTERFACE_H

#include <shape/BoundaryDiscretizer.hpp>

namespace imaging
{
  /** \ingroup shape 
      \brief Abstract class interface for discretizable shapes of dimension \em N.
      
      This class interface defines shapes of dimension \em N, which can be discretized. I.e. they provide the factory function boundary_discretizer(), which returns a boundary discretizer of dimension \em N for the given shape. A template function which integrates an image with floating point values along the boundary of an <em>N</em>-dimensional shape could generically be implemented as follows:
  \code
  template <size_t N>
  float_t boundary_area(const DiscretizableShapeInterface<N> & shape, const Image<N + 1, float_t> & image)
  {
    // sample 100 boundary points
    std::auto_ptr< BoundaryDiscretizer<N> > boundary_discretizer = shape.boundary_discretizer(100);
    
    return boundary_discretizer->integrate(image);
  }
  \endcode
  */
  template <size_t N>
  class DiscretizableShapeInterface
  {
  public:
    /** The spatial dimension of the shape class. \em N is 2 for planar shapes and 3 for volumes. */
    const static size_t SHAPE_DIMENSION = N;
    
    virtual ~DiscretizableShapeInterface() {};
    
    /** Returns a boundary discretizer for this shape. */ 
    virtual std::auto_ptr< BoundaryDiscretizer<SHAPE_DIMENSION> > boundary_discretizer(size_t n_points) const = 0;
  };
}


#endif
