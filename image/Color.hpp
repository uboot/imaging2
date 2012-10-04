// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_COLOR_H
#define IMAGE_COLOR_H

#include <core/imaging2.hpp>
#include <image/Image.hpp>
#include <image/GrayValue.hpp>

namespace imaging
{
  /** \ingroup image
      \brief Class for 24-bit color values.
      
      Objects of this class represent 24-bit color pixel values, which are internally stored as a triple of integers ranging from 0 to 255.
      
      \sa GrayValue
  */
  class Color : public ublas::fixed_vector<GrayValue, 3>
  {
  public:
    /** Corresponds to the RGB values (0, 0, 255). */
    static const Color BLUE;
    
    /** Corresponds to the RGB values (255, 0, 0). */
    static const Color RED;
    
    /** Corresponds to the RGB values (0, 255, 255). */
    static const Color CYAN;
    
    /** Corresponds to the RGB values (0, 255, 0). */
    static const Color GREEN;
    
    /** Corresponds to the RGB values (255, 255, 0). */
    static const Color YELLOW;
    
    /** Corresponds to the RGB values (255, 0, 255). */
    static const Color MAGENTA;
    
    /** Corresponds to the RGB values (0, 0, 0). */
    static const Color BLACK;
    
    /** Corresponds to the RGB values (255, 255, 255). */
    static const Color WHITE;

    Color() {}
    
    /** Copy constructor. */
    Color(const Color & color) : ublas::fixed_vector<GrayValue, 3>(color) {}
    
    /** Initializes object from \em value. This function sets each channel to \em value. */
    Color(GrayValue value) { (*this)(0) = value; (*this)(1) = value; (*this)(2) = value; }
    
    /** Initializes object from \em value. This function sets each channel to 255 * \em value. */
    Color(float_t value) { (*this)(0) = GrayValue(value); 
                            (*this)(1) = GrayValue(value);
                            (*this)(2) = GrayValue(value); }
                            
    /** Initializes object from red, green and blue values. */
    Color(unsigned char r, unsigned char g, unsigned char b)
      : ublas::fixed_vector<GrayValue, 3>(GrayValue(r), GrayValue(g), GrayValue(b)) {}
      
    /** Type cast to <em>float_t</em>. The color value is converted to GrayValue and then divided by 255.*/
    operator float_t() const { return float_t(GrayValue(*this)); }
  };
  
  /** \ingroup image
      \brief 2-dimensional color image.
      
      <tt>\#include <image/Color.hpp></tt>
      
      Abbreviated notation for 2-dimensional color images. The colors values are stored with 24 bit precision. i.e. each pixel is represented by 3 (red, green, blue) integer values 0-255.
  
      \sa Color
  */
  typedef Image<2, Color> ColorImage2d;
  
  template<>
  void ColorImage2d::read_image(const std::string & file_name);
  
  template<>
  void ColorImage2d::write_image(const std::string & file_name) const;
  
}

#endif
