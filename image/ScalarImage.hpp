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

#ifndef IMAGE_SCALARIMAGE_H
#define IMAGE_SCALARIMAGE_H

#include <image/ImageInterface.hpp>

namespace imaging
{

  /** \ingroup image
      \brief %Image with constant pixel values.
      
      This class models images with constant pixel values, i.e. every pixel has the same value. This value is set during construction. The size (in memory) of ScalarImage objects is independent of their size (as an image). 
  */
  template <std::size_t N, class DATA_t>
  class ScalarImage : public ImageInterface<N, DATA_t>
  {
  public:
    static const std::size_t dimension = N;
    typedef DATA_t data_t;

  private:
    ublas::fixed_vector<size_t, dimension> _size;
    data_t _value;

  public:
    /** Constructs a scalar image of size \em size. Each pixel of this image will evaluate to \em value. */
    ScalarImage(const ublas::fixed_vector<size_t, dimension> & size, data_t value) :
      _size(size), _value(value) {}
      
    /** Copy constructor. The complexity of this is function is constant (i.e. independent of the image size). */
    ScalarImage(const ScalarImage<dimension, data_t> & source) :
      _size(source._size), _value(source._value) {}
    
    const DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> &) const { return _value; }
    
    const ublas::fixed_vector<size_t, dimension> & size() const { return _size; }
  };


}

#endif
