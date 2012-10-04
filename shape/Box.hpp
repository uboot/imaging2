// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_BOX_H
#define SHAPE_BOX_H

#include <core/imaging2.hpp>
#include <core/utilities.hpp>

namespace imaging
{
  template <std::size_t N>
  class Box;
 
  /** \ingroup shape
      <tt>\#include <shape/Box.hpp></tt>
      
      Computes the intersection of two <em>N</em>-dimensional boxes and writes it to \em result.
  */
  template <std::size_t N>
  void compute_intersection_box(const Box<N> & box_a, const Box<N> & box_b, Box<N> & result);
    
   /** \ingroup shape
       \brief <em>N</em>-dimensional box (hyperrectangle).
       
       This class templates implements an <em>N</em>-dimensional hyperrectangle. For <em>N = 2</em> class objects are rectangles, for <em>N = 3</em> cuboids. The boxes are assumed to be aligned with the coordinate axes. I.e. each box is fully determined by two corner vertices.
  */
  template <std::size_t N>
  class Box
  {
    friend void compute_intersection_box<>(const Box<N> & box_a, const Box<N> & box_b, Box<N> & result);
  
    ublas::fixed_vector<float_t, N> _lower_corner;
    ublas::fixed_vector<float_t, N> _upper_corner;

  public:
    Box() : _lower_corner(), _upper_corner() {}
    
    /** Constructs a box from its two opposite corner vertices. */
    Box(const ublas::fixed_vector<float_t, N> & lower_corner, const ublas::fixed_vector<float_t, N> & upper_corner) : _lower_corner(lower_corner), _upper_corner(upper_corner)
    {}
    
    /** Constructs the bounding box of the image accessor \em accessor.*/
    template <class image_accessor_t>
    Box(const image_accessor_t & accessor)
    {
      for(size_t i = 0; i < N; ++i)
        _lower_corner = 0.0;
        
      _upper_corner = accessor.size();
    }
    
    /** Sets the box to the empty box. */
    void set_empty()
    {
      for(std::size_t i = 0; i < N; ++i)
      {
        _lower_corner(i) = 0.0;
        _upper_corner(i) = 0.0;
      }
    }
    
    /** Returns true if the box is empty. */
    bool is_empty() const
    {
      for(std::size_t i = 0; i < N; ++i)
        if(_lower_corner(i) >= _upper_corner(i))
          return true;
      
      return false;
    }
    
    /** Computes the volume of the box. */
    float_t compute_volume() const
    {
      float_t v = 1.0;
      
      if(is_empty())
        return 0.0;
        
      for(size_t i = 0; i < N; ++i)
        v *= _upper_corner(i) - _lower_corner(i);
        
      return v;
    }
    
    /** Determines if \em point lies within the box. If this is not the case the box is resized to the smalles box enclosing \em point and the original box. */
    void add_point(const ublas::fixed_vector<float_t, N> & point)
    {
      for(size_t i = 0; i < N; ++i)
      {
        if(point(i) > _upper_corner(i)) _upper_corner(i) = point(i);
        if(point(i) < _lower_corner(i)) _lower_corner(i) = point(i);
      }
    }
  };
  
  
  /** \ingroup shape
       <tt>\#include <shape/Box.hpp></tt>
      
       Determines if two <em>N</em>-dimensional boxes intersect.
  */
  template <std::size_t N>
  bool do_intersect(const Box<N> & box_a, const Box<N> & box_b)
  {
    if(box_a.is_empty() || box_b.is_empty())
      return false;
      
    Box<N> intersection_box;
    compute_intersection_box(box_a, box_b, intersection_box);
    
    if(intersection_box.is_empty())
      return false;
        
    return true;
  }
 
  template <std::size_t N>
  void compute_intersection_box(const Box<N> & box_a, const Box<N> & box_b, Box<N> & result)
  {
    for(size_t i = 0; i < N; ++i)
    {
      result._lower_corner(i) = max(box_a._lower_corner(i), box_b._lower_corner(i));
      result._upper_corner(i) = min(box_a._upper_corner(i), box_b._upper_corner(i));
    }
  }
  
  /** \ingroup shape
      <tt>\#include <shape/Box.hpp></tt>
      
      Computes the volume of the intersection of two <em>N</em>-dimensional boxes.
  */
  template <std::size_t N>
  float_t compute_intersection_volume(const Box<N> & box_a, const Box<N> & box_b)
  {
    if(! do_intersect(box_a, box_b))
      return 0.0;
      
    Box<N> intersection_box;
    
    compute_intersection_box(box_a, box_b, intersection_box);
    
    return intersection_box.compute_volume();
  }
}


#endif
