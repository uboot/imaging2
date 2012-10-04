// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
