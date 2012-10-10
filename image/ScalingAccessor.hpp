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

#ifndef IMAGE_SCALINGACCESSOR_H
#define IMAGE_SCALINGACCESSOR_H

#include <image/ImageAccessorInterface.hpp>

namespace imaging
{
  /** \ingroup image
      \brief Provides access to rescaled floating point images.
      
      This class allows you to access values of an image with floating point pixels after an affine transformation of the pixel values. Consider the following code:
  \code
  img::ScalingAccessor< img::Image<2, img::float_t> accessor(some_image, 1.0 / 128.0, -128.0);
  ublas::fixed_vector<2, img::size_t> lower_left(0, 0);
  float_t x = some_image[lower_left];
  float_t y = accessor[lower_left];
  accessor[lower_left] = 0.0;
  \endcode
      Then \f$x\f$ and \f$y\f$ are related in the following way:
      \f[
        y = \frac{1}{128} (x - 128) \quad\textrm{and}\quad x = 128 y + 128 \,.
      \f]
      I.e. \em accessor will scale an image with grayscale values in [0, 256[ to an image with values in [-1, 1[. If you assign a value in [-1, 1[ to the \em accessor it will implicitely convert it back to  [0, 255[ and assign this converted value to \em some_image. In the above example this means that <em>some_image[lower_left]</em> will be set to 128.
  */
  template <class float_accessor_t>
  class ScalingAccessor : public ImageAccessorInterface<float_accessor_t, float_t>
  {
  public:
    static const size_t dimension = ImageAccessorInterface<float_accessor_t, float_t>::dimension;
    typedef float_t data_t;

  private:
    /** \cond */
    class data_reference
    {
      ScalingAccessor<float_accessor_t> & _accessor_reference;
      const ublas::fixed_vector<size_t, dimension> & _index_reference;

    public:
      data_reference(ScalingAccessor<float_accessor_t> & accessor, const ublas::fixed_vector<size_t, dimension> & index) : _accessor_reference(accessor), _index_reference(index) {}

      float_t operator=(const float_t & input)
      { _accessor_reference._image_reference[_index_reference] = input / _accessor_reference._factor - _accessor_reference._level; return input; }

      operator float_t() const { return  _accessor_reference._factor * (_accessor_reference._image_reference[_index_reference] + _accessor_reference._level); }
    };
    /** \endcond */

    float_t _factor;
    float_t _level;

  public:
      /** Initialize the image accessor with an reference to \em image which must be an ImageInterface object storing floating point. This reference will be stored and accessed when the image accessor is utilized. The dimensions of the accessor are the same as those of \em image. Accessed pixel values in \em image will be first offset by \em level and then rescaled by \em factor.
      */
    ScalingAccessor(float_accessor_t & image, float_t factor = 1.0, float_t level = 0.0) :
    ImageAccessorInterface<float_accessor_t, float_t>(image), _factor(factor), _level(level) {}

    float_t operator[](const ublas::fixed_vector<size_t, dimension> & index) const { return _factor * (ImageAccessorInterface<float_accessor_t, float_t>::_image_reference[index] + _level); }
    data_reference operator[](const ublas::fixed_vector<size_t, dimension> & index) { return data_reference(*this, index); }
  };

}

#endif
