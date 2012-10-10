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

#ifndef IMAGE_GIO_H
#define IMAGE_GIO_H

#include <graphics/GraphicsInterface.hpp>
#include <image/Image.hpp>

namespace imaging
{
  /** \ingroup image 
      <tt>\#include <image/gio.hpp></tt>
      
      Draws \em image in the graphics output. The lower left corner of the image will be located at the origin (0, 0) and the upper right corner at (image.size()(0), image.size()(1)). In other words, this functions assumes the axes system of the graphics output to be scaled to pixel size. */
  GraphicsInterface & operator<<(GraphicsInterface & out, const ColorImage2d & image);
  
  template <std::size_t N, class DATA_t>
  GraphicsInterface & operator<<(GraphicsInterface & out, const Image<N, DATA_t> & image)
  {
    ColorImage2d colour_image(image.size());

    colour_image = image;
    
    out << colour_image;
    
    return out;
  }

}


#endif
