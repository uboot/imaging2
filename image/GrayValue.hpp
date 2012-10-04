// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_GRAYVALUE_H
#define IMAGE_GRAYVALUE_H

#include <core/imaging2.hpp>
#include <image/Image.hpp>

namespace imaging
{
  class Color;
  
  /** \ingroup image
      \brief Class for 8-bit grayscale values.
      
      Objects of this class represent 8-bit pixel values, which are internally stored as an integer from 0 to 255. If GrayValue objects are evalued as floating point objects (by calling <em>operator float_t()</em> or <em>GrayValue(float_t)</em> the value is scaled to the interval [0, 1].
      
      \sa Color
  */
  class GrayValue
  {
    unsigned char _value;
    
  public:
    GrayValue() {}
  
    /** Copy constructor. */
    GrayValue(const GrayValue & value) : _value(value._value) {}
    
    /** Construct a gray value from the colour value \em c_value. This function sets the gray value to the average of the red, green and blue values of \em c_value. */
    explicit GrayValue(const Color & c_value);

    /** Construct a gray value from \em value. */
    explicit GrayValue(unsigned char value) : _value(value) {}

    /** Construct a gray value from \em value. The input \em value is multiplied by 255 and then converted to an 8-bit value. */
    explicit GrayValue(float_t value) : _value((unsigned char)(value * 255.0)) {}

    /** Assignement operator. */
    GrayValue & operator=(GrayValue gray_value) { _value = gray_value._value; return *this; }

    /** Assignement operator. */
    GrayValue & operator=(unsigned char value) { *this = GrayValue(value); return *this; }

    /** Assignement operator. */
    GrayValue & operator=(float_t value) { *this = GrayValue(value); return *this; }
    
    /** Type cast to <em>unsigned char</em>. */
    operator unsigned char() const { return (unsigned char)(_value); }
    
    /** Type cast to <em>float_t</em>. The grayscale value is converted to floating point precision and divided by 255.*/
    operator float_t() const { return float_t(_value) / 255.0; }
  };
  
  /** \ingroup image
      \brief 2-dimensional grayscale image.
      
      <tt>\#include <image/GrayValue.hpp></tt>
      
      Abbreviated notation for 2-dimensional grayscale images. The gray values are stored with 8 bit precision. i.e. each pixel is represented by an integer value 0-255.
  
      \sa GrayValue
  */
  typedef Image<2, GrayValue> GrayImage2d;
  
  template<>
  void GrayImage2d::read_image(const std::string & file_name);
  
  template<>
  void GrayImage2d::write_image(const std::string & file_name) const;
  
  /** \ingroup image
      \brief 3-dimensional grayscale image.
      
      <tt>\#include <image/GrayValue.hpp></tt>
      
      Abbreviated notation for 3-dimensional grayscale images. The gray values are stored with 8 bit precision. i.e. each pixel is represented by an integer value 0-255.
  
      \sa GrayValue
  */
  typedef Image<3, GrayValue> GrayImage3d;
  
  template<>
  void GrayImage3d::read_image(const std::string & file_name);
  
  template<>
  void GrayImage3d::write_image(const std::string & file_name) const;
}

#endif
