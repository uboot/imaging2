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

#ifndef SHAPE_SHAPEINTERFACE_H
#define SHAPE_SHAPEINTERFACE_H

#include <core/imaging2.hpp>
#include <shape/BoundaryDiscretizer.hpp>

namespace imaging
{
  /** \ingroup shape 
      \brief Abstract class interface for shapes on a shape parameter manifold.
      
      This class interface defines shapes, which are determined by parameters on a Riemannian manifold. I.e. every shape implements the Riemannian exponential which maps vectors to a shape in its neighborhood. The Riemannian logarithm maps shapes in this neighborhood back to its Riemannian coordinates with respect to the reference shape. If we denote the Riemannian manifold as \f$M\f$ then every shape \f$p\f$ has to implement
      \f[
      \textrm{Exp}_p: T_pM \rightarrow M
      \f]
      and
      \f[
      \textrm{Log}_p: M \rightarrow T_pM
      \f]
      
      The following example shows how this enables computations with shapes without the knowledge of details about the shape implementation. It computes the Riemannian logarithm with respect to the first shape of each of the remaining shapes. Then the mean of the logarithms is computed and mapped back to a shape. The function is fully generic, i.e. it does not know what the actual shapes it processes are.
      
  \code
  template <class shape_t>
  void compute_mean(const std::vector<shape_t> & shapes, shape_t & mean_shape)
  {
    if(shapes.size() == 0)
      return;
      
    std::vector< ublas::vector<float_t> > logarithms(shapes.size() - 1);
    ublas::vector<float_t> logarithm_of_mean;
      
    shape_t & reference_shape = shapes[0]; // create an alias for the first shape
                                           // we call it reference_shape
                                           
    logarithm_of_mean = ublas::scalar_vector(reference_shape.dimension(), 0.0);
    // set the vector which will be the logarithm of the mean shape to zero
    
    ublas::vector<float_t> logarithm(reference_shape.dimension());
    // will hold the logarithm of each shape in the loop below
    
    for(size_t i = 1; i < shapes.size(); ++i)
    {
      reference_shape.logarithm(shapes[i], logarithm); // computes the logarithm
      logarithm_of_mean += logarithm; // adds it to the scaled mean
    }
    
    logarithm_of_mean /= float_t(shapes.size()); // scales the mean correctly
    
    reference_shape.exponential(logarithm_of_mean, mean_shape);
    // computes the shape determined by logarithm_of_mean
  }
  \endcode
  */
  class ShapeInterface
  {
  public:
    virtual ~ShapeInterface() {};
    
    /** Computes the Riemannian exponential of \em vector with respect to \em *this and stores the result in \em shape. In other words, \em shape is the projection of the tangent vector \em vector in the tangent space at the \em *this shape on the shape manifold. */
    virtual void exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const = 0;
    
    /** Computes the Riemannian logarithm of \em shape with respect to \em *this and stores the result in \em vector. In other words, \em vector is the projection of \em shape on the tangent space at the \em *this shape on the shape manifold. */
    virtual void logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const = 0;
    
    /** Returns the parametric dimension of the shape object. This is the same as the dimension of the shape manifold the shape belongs to. It should not be confused with the spatial dimension \em N. */
    virtual size_t dimension() const = 0;
  };
}


#endif
