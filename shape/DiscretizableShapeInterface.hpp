// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
