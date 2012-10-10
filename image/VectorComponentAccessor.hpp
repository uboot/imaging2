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

#ifndef VECTORCOMPONENTACCESSOR_H
#define VECTORCOMPONENTACCESSOR_H

#include <core/imaging2.hpp>

namespace imaging
{

  /** \ingroup image
      \brief Provides access to different channels of vector valued image data.
      
      Use this class to access components of vector valued pixels in images. By <em>vector valued</em> we mean that the pixel data type of the input image data consists of different components (all of the same type) which can be accessed by <em>%operator()</em>. The following code example copies the green component of \em colour_image to \em green_channel. It then processes the green_channel and assigns the possibly modified channel back to the original image:
  \code
  ColorImage2d colour_image("some_image.png");
  GrayImage2d green_channel(colour_image.size());
  
  VectorComponentAccessor<ColorImage2d> accessor(colour_image, 2);
  
  green_channel = accessor;
  
  // process green_channel...
  
  accessor = green_channel;
  \endcode
  */
  template <class vector_image_t>
  class VectorComponentAccessor : public ImageAccessorInterface<vector_image_t, typename vector_image_t::data_t::data_t>
  {
  public:
    typedef typename vector_image_t::data_t::data_t data_t;
    static const size_t dimension =
      ImageAccessorInterface<vector_image_t, data_t>::dimension;
  private:
    /** \cond */
    class data_reference
    {
      VectorComponentAccessor<vector_image_t> & _accessor_reference;
      const ublas::fixed_vector<size_t, dimension> & _index_reference;

    public:
      data_reference(VectorComponentAccessor<vector_image_t> & accessor, const ublas::fixed_vector<size_t, dimension> & index) : _accessor_reference(accessor), _index_reference(index) {}

      float_t operator=(const float_t & input)
      { _accessor_reference._image_reference[_index_reference](_accessor_reference._component) = input; return input; }

      operator float_t() const { return _accessor_reference._image_reference[_index_reference](_accessor_reference._component); }
    };
    /** \endcond */

    size_t _component;

  public:
    /** Initialize the image accessor with an reference to \em image which must be an ImageInterface object with vector valued pixel data. The accessor will then behave as a scalar valued image consisting of the values of the <em>component</em>-th channel of \em image. The reference will be stored and accessed when the image accessor is utilized. The dimensions of the accessor are the same as those of \em image.
        
        Note that \em component must be within the range of the components of the pixel data of \em image. E.g. \em component must not exceed 2 in case of RGB images.
    */
    VectorComponentAccessor(vector_image_t & image, size_t component) : ImageAccessorInterface<vector_image_t, data_t>(image), _component(component) {}

    data_t operator[](const ublas::fixed_vector<size_t, dimension> & index) const { return ImageAccessorInterface<vector_image_t, data_t>::_image_reference[index](_component); }

    data_reference operator[](const ublas::fixed_vector<size_t, dimension> & index) { return data_reference(*this, index); }
  };
}

#endif
