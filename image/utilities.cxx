#include <image/utilities.hpp>

#include <core/utilities.hpp>
#include <Magick++.h>


namespace imaging
{
  template <>
  void edge(Image<2, float_t> & image, float_t radius)
  {
    Magick::Image magick_image(image.size()(1), image.size()(0), "K", Magick::DoublePixel, &image[ublas::fixed_vector<size_t, 2>(0, 0)]);
    
    magick_image.edge(radius);
    
    magick_image.write(0, 0, image.size()(1), image.size()(0), "K", Magick::DoublePixel, &image[ublas::fixed_vector<size_t, 2>(0, 0)]);
  }
  
  template <>
  void blur(Image<2, float_t> & image, float_t width, float_t sigma)
  {
    if(sigma == 0.0) return;
    
    Magick::Image magick_image(image.size()(1), image.size()(0), "K", Magick::DoublePixel, &image[ublas::fixed_vector<size_t, 2>(0, 0)]);
    
    magick_image.gaussianBlur(width, sigma);
    
    magick_image.write(0, 0, image.size()(1), image.size()(0), "K", Magick::DoublePixel, &image[ublas::fixed_vector<size_t, 2>(0, 0)]);
  }
  
  template <>
  void absolute_gradient(Image<2, float_t> & image)
  {
    Image<2, float_t> temp_image(image.size());
    
    ublas::fixed_vector<size_t, 2> h(1, 0), v(0, 1);
    
    for(size_t i = 1; i < image.size()(0) - 1; ++i)
      for(size_t j = 1; j < image.size()(1) - 1; ++j)
      {
        ublas::fixed_vector<size_t, 2> index(i, j);
        temp_image[index] = sqrt(square( (image[index - h] - 2 * image[index] + image[index + h]) / 2.0 ) + 
                                 square( (image[index - v] - 2 * image[index] + image[index + v]) / 2.0 ) );
      }
      
          
    for(size_t i = 1; i < image.size()(0) - 1; ++i)
    {
      ublas::fixed_vector<size_t, 2> index(i, 0);
      temp_image[index] = temp_image[index + v];
                               
      index.assign(i, image.size()(1) - 1);
      temp_image[index] = temp_image[index - v];
    }
      
       
    for(size_t j = 1; j < image.size()(1) - 1; ++j)
    {
      ublas::fixed_vector<size_t, 2> index(0, j);
      temp_image[index] = temp_image[index + h];
                               
      index.assign(image.size()(0) - 1, j);
      temp_image[index] = temp_image[index - h];
    }
    
    ublas::fixed_vector<size_t, 2> index(0, 0);
    temp_image[index] = ( temp_image[index + h] + temp_image[index + v] ) / 2.0;
                             
    index.assign(0, image.size()(1) - 1);
    temp_image[index] = ( temp_image[index + h] + temp_image[index - v] ) / 2.0;
                             
    index.assign(image.size()(0) - 1, 0);
    temp_image[index] = ( temp_image[index - h] + temp_image[index + v] ) / 2.0;
                             
    index.assign(image.size()(0) - 1, image.size()(1) - 1);
    temp_image[index] = ( temp_image[index - h] + temp_image[index - v] ) / 2.0;
    
      
    image = temp_image;
  }
  
  void draw_level_set(const Image<2, float_t> & level_set_function, const float_t & level, const Color & color, ColorImage2d & image)
  {
    if(level_set_function.size() != image.size())
      throw Exception("Exception: Dimensions of level set function and image do not agree filter_plugin::geodesic_active_contours::draw_level_set().");
      
    for(size_t i = 1; i < level_set_function.size()(0) - 1; ++i)
      for(size_t j = 1; j < level_set_function.size()(1) - 1; ++j)
      {
        float_t center = level_set_function[ublas::fixed_vector<size_t, 2>(i, j)] - level;
        float_t right = level_set_function[ublas::fixed_vector<size_t, 2>(i + 1, j)] - level;
        float_t left = level_set_function[ublas::fixed_vector<size_t, 2>(i - 1, j)] - level;
        float_t top = level_set_function[ublas::fixed_vector<size_t, 2>(i, j + 1)] - level;
        float_t bottom = level_set_function[ublas::fixed_vector<size_t, 2>(i, j - 1)] - level;
        
        if(center * right < 0.0 || center * left < 0.0 || center * top < 0.0 || center * bottom < 0.0)
          image[ublas::fixed_vector<size_t, 2>(i, j)] = color;
      }
  }
  
  void image2matrix(const Image<2, float_t> & image, ublas::matrix<float_t> & matrix)
  {
    matrix.resize(image.size()(0), image.size()(1));
    
    for(size_t i = 0; i < image.size()(0); ++i)
      for(size_t j = 0; j < image.size()(1); ++j)
        matrix(i, j) = image[ublas::fixed_vector<float_t, 2>(i, j)];
  }
}

