// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
