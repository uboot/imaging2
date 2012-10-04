// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_INTERPOLATIONADAPTORINTERFACE_H
#define IMAGE_INTERPOLATIONADAPTORINTERFACE_H

// TODO: fix handling of images of size 1
namespace imaging
{
  namespace interpolation_adaptor_impl
  {
    template <size_t N>
    bool offset_position(const ublas::fixed_vector<size_t, N> & size, const ublas::fixed_vector<float_t, N> & position, ublas::fixed_vector<float_t, N> & out)
    {      
      for(size_t i = 0; i < N; ++i)
      {
        if(position(i) < 0.0 || position(i) >= float_t(size(i)) || size(i) < 2)
        {
          return false;
        }
        
        if(position(i) <= 0.5)
          out(i) = 0.0;
        else if(position(i) >= float_t(size(i)) - 0.5)
          out(i) = float_t(size(i) - 1);
        else 
          out(i) = position(i) - 0.5;
      }
      
      return true;
    }
    
    template <size_t N>
    void compute_local_coordinates(const ublas::fixed_vector<float_t, N> & index, const ublas::fixed_vector<size_t, N> & size, ublas::fixed_vector<size_t, N> & base_coordinates, ublas::fixed_vector<float_t, N> & truncated_coordinates)
    { 
      for(size_t i = 0; i < N; ++i)
      {
        base_coordinates(i) = size_t(floor(index(i)));
        truncated_coordinates(i) = index(i) - float_t(base_coordinates(i));
        
        if(base_coordinates(i) == size(i) - 1)
        {
          base_coordinates(i)--;
          truncated_coordinates(i) = 1.0;
        }
      }
    }
  }
  
  template <class image_t>
  class PolynomialInterpolationAdaptor;
  
  template <class image_t>
  class LinearInterpolationAdaptor;
  
  template <class image_t>
  class MeanInterpolationAdaptor;

  /** \ingroup image
      \brief Abstract base class of image interpolation function which provide pixel values in floating point coordinates.
      
      This class is the base class of all image interpolation classes.
      Such classes are initialized by ImageInterface objects and provide a pixel accessor [] for real valued pixel coordinates.
      Depending on the implementing class this operators returns an interpolated pixel value.
  */
  template <class image_t>
  class InterpolationAdaptorInterface
  {
    friend class PolynomialInterpolationAdaptor<image_t>;
    friend class LinearInterpolationAdaptor<image_t>;
    friend class MeanInterpolationAdaptor<image_t>;
    
    const image_t & _image_reference;
  
    InterpolationAdaptorInterface(const image_t & image_reference) : _image_reference(image_reference)
    {}

  public:
    /** The data type which is returned by the interpolation adaptor. Common choices are \em float_t for grayscale images stored at floating point precision or \em char for 8-bit data. */
    typedef typename image_t::data_t data_t;
    
    /** The dimension of the image. It equals 2 for planar images and 3 for voxel-data. */
    static const std::size_t dimension;
    
    /** Returns the size of the image the interpolation adaptor is attached to. */
    const ublas::fixed_vector<size_t, image_t::dimension> & size() const
    {
      return _image_reference.size();
    }
    
    /** Returns the interpolation at the position \em index of the surrounding pixel values. */
    data_t operator[](const ublas::fixed_vector<float_t, image_t::dimension> & index);
  };
  
  template <class image_t>
  const std::size_t InterpolationAdaptorInterface<image_t>::dimension = image_t::dimension;
}

#endif
