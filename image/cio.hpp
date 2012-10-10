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

#ifndef IMAGE_CIO_H
#define IMAGE_CIO_H


#include <image/ImageInterface.hpp>
#include <image/GrayValue.hpp>
#include <image/Color.hpp>


namespace imaging
{
  /** \ingroup image
      <tt>\#include <image/cio.hpp></tt>
  */
  template <size_t N, class DATA_t>
  std::ostream & operator<<(std::ostream & out, const ImageInterface<N, DATA_t> & image)
  {
    out << "Image: " << image.size()(0);
    
    for(size_t i = 1; i < ImageInterface<N, DATA_t>::dimension; ++i)
      out << " x " << image.size()(i);
      
    out << std::endl;
      
    return out;
  }
 
  /** \ingroup image
      <tt>\#include <image/cio.hpp></tt>
  */
  std::ostream & operator<<(std::ostream & out, const GrayValue value);
  
  /** \ingroup image
      <tt>\#include <image/cio.hpp></tt>
  */
  std::ostream & operator<<(std::ostream & out, const Color & value);

  /** \ingroup image
      <tt>\#include <image/cio.hpp></tt>
  */
  std::istream & operator>>(std::istream & in, GrayValue & value);
  
  /** \ingroup image
      <tt>\#include <image/cio.hpp></tt>
  */
  std::istream & operator>>(std::istream & in, Color & value);
}



#endif
