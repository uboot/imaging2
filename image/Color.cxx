#include <image/Color.hpp>

#include <stdio.h>
#include <Magick++.h>

#include <core/MessageInterface.hpp>


namespace imaging
{
  const Color Color::BLUE(0, 0, 255);
  const Color Color::RED(255, 0, 0);
  const Color Color::CYAN(0, 255, 255);
  const Color Color::GREEN(0, 255, 0);
  const Color Color::YELLOW(255, 255, 0);
  const Color Color::BLACK(0, 0, 0);
  const Color Color::WHITE(255, 255, 255);
  const Color Color::MAGENTA(255, 0, 255);
  
  template<>
  void ColorImage2d::read_image(const std::string & file_name)
  {
    Magick::Image magick_image;
    
    try
    { 
      magick_image.read(file_name); 
    } 
    catch(Magick::Exception &error_) 
    { 
      throw FileIoException
      ("FileIoException: File " + file_name + " could not be opened for reading in GrayImage2d::read_image().");
    } 
    
    ColorImage2d color_image(ublas::fixed_vector<size_t, 2>(magick_image.size().height(), magick_image.size().width()));
    magick_image.write(0, 0, magick_image.size().width(), magick_image.size().height(), "RGB", Magick::CharPixel, &color_image[ublas::fixed_vector<size_t, 2>(0, 0)]);
    
    this->resize(ublas::fixed_vector<size_t, 2>(color_image.size()(1), color_image.size()(0)));
    
    for(size_t i = 0; i < color_image.size()(0); ++i)
      for(size_t j = 0; j < color_image.size()(1); ++j)
        (*this)[ublas::fixed_vector<size_t, 2>(j, color_image.size()(0) - i - 1)] = color_image[ublas::fixed_vector<size_t, 2>(i, j)];
  }
  
  template<>
  void ColorImage2d::write_image(const std::string & file_name) const
  {
    ColorImage2d rotated_image(ublas::fixed_vector<size_t, 2>(this->size()(1), this->size()(0)));
    
    for(size_t i = 0; i < rotated_image.size()(0); ++i)
      for(size_t j = 0; j < rotated_image.size()(1); ++j)
        rotated_image[ublas::fixed_vector<size_t, 2>(i, j)] = (*this)[ublas::fixed_vector<size_t, 2>(j, this->size()(1) - i - 1)];
    
    Magick::Image magick_image(rotated_image.size()(1), rotated_image.size()(0), "RGB", Magick::CharPixel, &rotated_image[ublas::fixed_vector<size_t, 2>(0, 0)]);
    
    magick_image.write(file_name);
  }
}

