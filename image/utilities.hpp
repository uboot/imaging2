// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_UTILITIES_H
#define IMAGE_UTILITIES_H

#include <image/Image.hpp>
#include <image/Color.hpp>

namespace imaging
{
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Computes a vector field \em out such that the divergence of \em out equals the scalar input \em image. If we denote \image as \f$f\f$ and \em out as \f$\vec F\f$ then this means that
      \f[
        \nabla \cdot \vec F = f\,.
      \f]
      This is done by setting
      \f[
        \vec F_n(x_1, \ldots, x_n, \ldots, x_N) = \frac1N \int_0^{x_n}f(x_1, \ldots, \xi, \ldots, x_N) d\xi\,,
      \f]
      where \f$N\f$ denotes the dimension of \em image. 
      
      The dimension of \em image must equal the number of components of each pixel in \em out. 
  */
  template <class float_accessor_t, class vector_accessor_t>
  void compute_divergence_field(const float_accessor_t & image, vector_accessor_t & out)
  {
    const size_t dimension = float_accessor_t::dimension;

    if(image.size() != out.size())
      throw Exception("Exception: Dimensions in compute_divergence_field() do not agree.");

    for(size_t i = 0; i < dimension; ++i)
    {
      ublas::fixed_vector<size_t, dimension> index;
      index.assign(0);

      if(image.size()(i) == 0)
        break;

      do
      {
        out[index](i) = 1.0/float_t(dimension) * image[index];
      }
      while(increment_index(image.size(), i, index));

      ublas::fixed_vector<size_t, dimension> offset;
      offset.assign(0);

      for(size_t j = 1; j < image.size()(i); ++j)
      {
        index(i) = j;
        do
        {
          offset(i) = 1;
          out[index](i) = out[index - offset](i) + 1.0/float_t(dimension) * image[index];
          offset(i) = 0;
        }
        while(increment_index(image.size(), i, index));
      }
    }
  }
  
  /** \ingroup image 
      <tt>\#include <image/utilities.hpp></tt>
      
      Blurs \em image by convoluting it with a Gaussian kernel with standard deviation \em sigma (in pixels). The parameter \em width determines the radius of the kernel in pixels. This function is currently implemented for 2-dimensional images only (using ImageMagick).
  */
  template <class image_t>
  void blur(image_t & image, float_t width, float_t sigma)
  {}
  
  template <>
  void blur(Image<2, float_t> & image, float_t width, float_t sigma);
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Calls the \em edge function of ImageMagick. No idea what it does. Implemented for 2-dimensional images only.
  */
  template <class image_t>
  void edge(image_t & image, float_t radius)
  {}
  
  template <>
  void edge(Image<2, float_t> & image, float_t radius);
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
        Computes the absolute value of the gradient of \em image in the following way:
          - the gradient of the inner pixels of \em image is computed using symmetric finite differences,
          - the inner pixels of \em image are set to the absolute value of the gradients computed in the previous step,
          - the values of the boundary pixels (excluding the pixels in the corners) \em image are set to the values of the neighboring inner pixels,
          - the values of the boundary pixels are set to the mean of the values in the neighboring pixels.
        This is implemented for 2-dimensional images only, but the above procedure generalizes to higher dimensions (somebody must implement it, though).
  */
  template <class image_t>
  void absolute_gradient(image_t & image)
  {}
  
  template <>
  void absolute_gradient(Image<2, float_t> & image);
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Computes the maximal pixel value of \em image. The class \em image_t must implement ImageInterface and \em image_t::data_t must be comparable.
  */
  template <class image_t>
  typename image_t::data_t max(const image_t & image)
  {
    if(image.empty()) 
      throw Exception("Exception: passed empty image to max(const ImageInterface &).");
    
    ublas::fixed_vector<size_t, image_t::dimension> index;
    index.assign(0);
    typename image_t::data_t max_value = image[index];
    do
    {
      typename image_t::data_t temp = image[index];
      if(max_value < temp) max_value = temp;
    }
    while(increment_index(image.size(), index));
    
    return max_value;
  }
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Computes the minimal pixel value of \em image. The class \em image_t must implement ImageInterface and \em image_t::data_t must be comparable.
  */
  template <class image_t>
  typename image_t::data_t min(const image_t & image)
  {
    if(image.empty()) 
      throw Exception("Exception: passed empty image to max(const ImageInterface &).");
    
    ublas::fixed_vector<size_t, image_t::dimension> index;
    index.assign(0);
    typename image_t::data_t min_value = image[index];
    do
    {
      typename image_t::data_t temp = image[index];
      if(min_value > temp) min_value = temp;
    }
    while(increment_index(image.size(), index));
    
    return min_value;
  }
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Returns true if the pixel defined by \em index is within the bounds of \em image and false otherwise.
  */
  template <class image_t>
  bool is_pixel(const image_t & image, const ublas::fixed_vector<size_t, image_t::dimension> & index)
  {
    for(size_t i = 0; i < image_t::dimension; ++i)
      if(index(i) >= image.size()(i))
        return false;
    
    return true;
  }
  
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Draws the level line of \em level_set_function defined by \em level onto \em image using \em color.
      Implemented for 2-dimensional images only.
  */
  void draw_level_set(const Image<2, float_t> & level_set_function, const float_t & level, const Color & color, ColorImage2d & image);
  
  /** \ingroup image
      <tt>\#include <image/utilities.hpp></tt>
      
      Converts the 2-dimensional \em image to \em matrix.
  */
  void image2matrix(const Image<2, float_t> & image, ublas::matrix<float_t> & matrix);
}

#endif
