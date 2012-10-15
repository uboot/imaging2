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

#ifndef SHAPE_BSPLINESHAPE_H
#define SHAPE_BSPLINESHAPE_H

#include <shape/ShapeInterface.hpp>
#include <shape/DiscretizableShapeInterface.hpp>

#include <spline/PeriodicBspline.hpp>

namespace imaging
{
  /** \ingroup shape
      \brief Class for shapes determined by boundary spline curves in the plane.
      
      This class implements shapes in the plane which are determined by a B-spline boundary curve.
  */
  class BsplineShape : public ShapeInterface, public DiscretizableShapeInterface<2>
  {
    PeriodicBspline< ublas::fixed_vector<float_t, 2> > _curve;
    
    class Discretizer;
    
  public:
    const static size_t SHAPE_DIMENSION = 2;

    BsplineShape() : _curve() {}
    
    /** Copy constructor. */
    BsplineShape(const BsplineShape & source) : _curve(source._curve) {}
    
    /** Constructs a B-spline shape from the a B-spline \em curve. */
    BsplineShape(const PeriodicBspline< ublas::fixed_vector<float_t, 2> > & curve) : _curve(curve) {}

    /** Copy assignement. */
    const BsplineShape & operator=(const BsplineShape & source)
    { _curve = source._curve; return *this; }
    
    virtual std::auto_ptr< BoundaryDiscretizer<2> > boundary_discretizer(size_t n_points) const;

    void exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const;

    void logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const;

    size_t dimension() const { return 2 * _curve.n_coefficients(); }
    
    /** Returns the boundary curve of the shape. */
    const PeriodicBspline< ublas::fixed_vector<float_t, 2> > & curve() const { return _curve; }
  };
  
  
  /** \cond */
  class BsplineShape::Discretizer : public BoundaryDiscretizer<2>
  {
    float_t _step_size;
    const BsplineShape & _spline_curve;

  public:
  
    Discretizer(const BsplineShape  & spline_curve, size_t n_points);
    
    void evaluate(size_t i, ublas::fixed_vector<float_t, 2> & point, ublas::fixed_vector<float_t, 2> & normal, float_t & curvature) const;
  };
  /** \endcond */
  
}


#endif
