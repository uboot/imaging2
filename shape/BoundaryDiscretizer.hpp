// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_BOUNDARYDISCRETIZER_H
#define SHAPE_BOUNDARYDISCRETIZER_H

#include <shape/Box.hpp>

namespace imaging
{
  /** \cond */
  namespace boundary_discretizer_impl
  {
    float_t compute_curve_curvature(ublas::fixed_vector<float_t, 2> & first_derivative, ublas::fixed_vector<float_t, 2> & second_derivative);
    
    template <class const_vector_accessor_t>
    void compute_zero_extension(const const_vector_accessor_t  & vector_field, const ublas::fixed_vector<float_t, const_vector_accessor_t::dimension> & point, ublas::fixed_vector<float_t, const_vector_accessor_t::dimension> & value)
    {
      const size_t N = const_vector_accessor_t::dimension;
      bool is_valid = true;
      bool no_valid_index = true;
      ublas::fixed_vector<size_t, N> valid_index;
      value = ublas::scalar_vector<float_t>(N, 0.0);
      
      // check if a coordinate is negative
      for(size_t i = 0; i < N; ++i)
        if(point(i) < 0.0)
          return;
        
      valid_index = ublas::fixed_vector<size_t, N>(point);
      int valid_coordinate;

      for(size_t i = 0; i < N; ++i)
        if( ! ( valid_index(i) < vector_field.size()(i) ) )
        {
          if(! is_valid)
            return;
          else  
          {
            valid_index(i) = vector_field.size()(i) - 1;
            is_valid = false;
          }
          valid_coordinate = i;
          no_valid_index = false;
        }
      
      if(no_valid_index)
      {
        value = vector_field[ublas::fixed_vector<size_t, N>(point)];
        return;
      }
      
      value(valid_coordinate) = vector_field[valid_index](valid_coordinate);
    }
  }
  
  /** \endcond */
  /** \ingroup shape 
      \brief Abstract class interface for discretizations of shape boundaries. 
      
      Objects derived from this class can be evaluated in a finite number of sample points on the boundary of a shape. In addition to the coordinates of the sample point they must provide the outer normal in that point. The length of the normal must be such that the lengthes of all normals sum up to the area of the shape boundary. In other words, the length of the normal corresponds to the area of the infinitesimal boundary element around the normal.
  */
  template <size_t N>
  class BoundaryDiscretizer
  {
  protected:
    /** Number of discretization points. */
    size_t _n_points;
    
  public:
    /** The spatial dimension of the shape boundary. \em N is 2 for planar shapes and 3 for volumes. */
    static const size_t SHAPE_DIMENSION = N;
  
    /** Construct a discretization with \em n_points discretization points. */ 
    BoundaryDiscretizer(size_t n_points) : _n_points(n_points) {}
    
    virtual ~BoundaryDiscretizer() {}
    
    /** Returns the number of discretization points. */
    size_t n_points() const { return _n_points; }
    
    /** Sets \em point to the coordinates of the <em>i</em>-th disretization point and stores the boundary normal in this point in \em normal. The normal must be scaled to the area of the infinitesimal boundary element around the point. This function must be implemented in classes derived from BoundaryDiscretizer. For the actual evaluation of a discretization point the user should use operator(). */ 
    virtual void evaluate(size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & point, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const = 0;
    
    /** Evaluates the coordinates of the <em>i</em>-th discretization point, sets \em normal to the boundary normal and \em curvature to the curvature in this point. The normal is scaled to the area of the infinitesimal boundary element around the point. */
    ublas::fixed_vector<float_t, SHAPE_DIMENSION> operator()(size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const
    {
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> point;
      evaluate(i, point, normal, curvature);
      return point;
    }
    
    /** Evaluates the coordinates of the <em>i</em>-th discretization point and sets \em normal to the boundary normal in this point. The normal is scaled to the area of the infinitesimal boundary element around the point. */
    ublas::fixed_vector<float_t, SHAPE_DIMENSION> operator()(size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal) const
    {
      float_t curvature;
      return operator()(i, normal, curvature);
    }
    
    /** Evaluates the coordinates of the <em>i</em>-th discretization point. */
    ublas::fixed_vector<float_t, SHAPE_DIMENSION> operator()(size_t i) const
    {
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> normal;
      return operator()(i, normal);
    }
      
    /** Integrate an image over a boundary discretization. */
    template <class const_accessort_t>
    typename const_accessort_t::data_t integrate(const const_accessort_t & image) const
    {
      typename const_accessort_t::data_t result = 0.0;
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> point, normal;

      for(size_t i = 0; i < n_points(); ++i)
      {
        point = (*this)(i, normal);
        
        // explicit cast of floating point values in "vertex" to pixel position!
        ublas::fixed_vector<size_t, const_accessort_t::dimension> pixel_position = ublas::fixed_vector<size_t, const_accessort_t::dimension>(point);
        bool is_valid_pixel = true;
        
        for(size_t i = 0; i < const_accessort_t::dimension; ++i)
        {
          if( ! ( pixel_position(i) < image.size()(i) ) )
          {
            is_valid_pixel = false;
            break;
          }
        }
        
        if(is_valid_pixel)
          result += image[pixel_position] * norm_2(normal);
      }

      return result;
    }
    
    /** Compute the <em>(n-1)</em>-dimensional boundary area of a boundary discretization. */
    float_t compute_boundary_area() const
    {
      float_t result = 0.0;
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> point, normal;

      for(size_t i = 0; i < n_points(); ++i)
      {
        point = (*this)(i, normal);
        
        result += norm_2(normal);
      }

      return result;
    }
    
    
    /** Integrate a vector field over a boundary discretization. */
    template <class const_vector_accessor_t>
    float_t integrate_vector_field
    (const const_vector_accessor_t & vector_field) const
    {
      float_t result = 0.0;
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> point, normal;

      for(size_t i = 0; i < n_points(); ++i)
      {
        point = (*this)(i, normal);
        
        // explicit cast of floating point values in "vertex" to pixel position!
        ublas::fixed_vector<size_t, const_vector_accessor_t::dimension> pixel_position = ublas::fixed_vector<size_t, const_vector_accessor_t::dimension>(point);
        bool is_valid_pixel = true;
        
        for(size_t i = 0; i < const_vector_accessor_t::dimension; ++i)
        {
          if( ! ( pixel_position(i) < vector_field.size()(i) ) )
          {
            is_valid_pixel = false;
            break;
          }
        }
        
        if(is_valid_pixel)
          // Note: unit normal * tangential speed = normal
          result += inner_prod(vector_field[pixel_position], normal);
        else
        {
          ublas::fixed_vector<float_t, SHAPE_DIMENSION> value;
          boundary_discretizer_impl::compute_zero_extension(vector_field, point, value);
          result += inner_prod(value, normal);
        }
      }

      return result;
    }
    
            
    
    /** Compute the bounding box of a boundary discretization. */
    Box<SHAPE_DIMENSION> compute_bounding_box() const
    {
      ublas::fixed_vector<float_t, SHAPE_DIMENSION> lower_corner, upper_corner, point;
      
      if(n_points() == 0)
        return Box<SHAPE_DIMENSION>();
        
      lower_corner = (*this)(0);
      upper_corner = (*this)(0);
          
      for(size_t i = 1; i < n_points(); ++i)
      {
        point = (*this)(i);
        for(size_t j = 0; j < SHAPE_DIMENSION; ++j)
        {
          if(point(j) < lower_corner(j))
            lower_corner(j) = point(j);
          
          if(point(j) > upper_corner(j))
            upper_corner(j) = point(j);
        }
      }
      
      return Box<SHAPE_DIMENSION>(lower_corner, upper_corner);
    }
  };
  
  /** Compute the area of intersection of two (2-dimensional) discretized shapes. */
  float_t compute_intersection_volume(const BoundaryDiscretizer<2> & shape_a, const BoundaryDiscretizer<2> & shape_b);    
}


#endif
