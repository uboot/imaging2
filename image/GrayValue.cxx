#include <image/GrayValue.hpp>

#include <image/Color.hpp>
#include <core/MessageInterface.hpp>

#include <stdio.h>
#include <stdlib.h>


namespace imaging
{
  GrayValue::GrayValue(const Color & c_value)
  {
    unsigned int temp = 0;

    for(size_t i = 0; i < 3; ++i)
      temp += (unsigned char)(c_value(i));

    _value = temp/3;
  }
  
  template<>
  void GrayImage2d::read_image(const std::string & file_name)
  {
    ColorImage2d color_image(file_name);
      
    *this = color_image;
  }
  
  
  template<>
  void GrayImage2d::write_image(const std::string & file_name) const
  {
    ColorImage2d color_image(*this);
    
    color_image.write_image(file_name);
  }
  
  
  template<>
  void GrayImage3d::read_image(const std::string & file_name)
  {
    FILE * in_file; 
    long file_size;
    uint32_t size_buffer[3];
    ublas::fixed_vector<size_t, 3> image_size;
    size_t result;
    
    in_file = fopen (file_name.c_str() , "rb" );
    if (in_file == NULL) 
      throw FileIoException("FileIoException: could not open file '" + file_name +"' in GrayImage3d::read_image().");
    
    fseek(in_file , 0 , SEEK_END);
    file_size = ftell (in_file);
    rewind (in_file);
    
    if(file_size < sizeof(uint32_t) * 3)
      throw FileIoException("FileIoException: file '" + file_name +"' too short in GrayImage3d::read_image().");
      
    result = fread(size_buffer, 1, sizeof(uint32_t) * 3 , in_file); 
    
    if (result != sizeof(uint32_t) * 3)
      throw FileIoException("FileIoException: read error GrayImage3d::read_image().");
    
    for(size_t i = 0; i < 3; ++i)
      image_size(i) = size_buffer[i];
    
    resize(image_size);
    
    for(size_t z = 0; z < image_size(2); ++z)
      for(size_t y = 0; y < image_size(1); ++y)
        for(size_t x = 0; x < image_size(0); ++x)
        {
          int c = fgetc(in_file);
          if(c == EOF)
            throw FileIoException("FileIoException: file '" + file_name +"' too short in GrayImage3d::read_image().");
          (*this)[ublas::fixed_vector<size_t, 3>(x, y, z)] = (unsigned char)(c);
        }
        
    fclose(in_file);
  }
  
  template<>
  void GrayImage3d::write_image(const std::string & file_name) const
  {
    FILE * out_file; 
    uint32_t size_buffer[3];
    size_t result;
    
    for(size_t i = 0; i < 3; ++i)
      size_buffer[i] = size()(i);
    
    out_file = fopen (file_name.c_str() , "wb" );
    
    if (out_file == NULL) 
      throw FileIoException("FileIoException: could not open file '" + file_name +"' in GrayImage3d::write_image().");
      
    result = fwrite(size_buffer, 1, sizeof(uint32_t) * 3, out_file);
    
    if (result != sizeof(uint32_t) * 3)
      throw FileIoException("FileIoException: write error GrayImage3d::write_image().");
      
    for(size_t z = 0; z < size()(2); ++z)
      for(size_t y = 0; y < size()(1); ++y)
        for(size_t x = 0; x < size()(0); ++x)
        {
          int c = fputc((unsigned char)(*this)[ublas::fixed_vector<size_t, 3>(x, y, z)], out_file);
          if(c != (unsigned char)(*this)[ublas::fixed_vector<size_t, 3>(x, y, z)])
            throw FileIoException("FileIoException: write error GrayImage3d::write_image().");
        }  
    
    fclose(out_file);     
  }
  
}

