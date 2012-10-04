// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_POLYNOMIALINTERPOLATIONADAPTOR_H
#define IMAGE_POLYNOMIALINTERPOLATIONADAPTOR_H

#include <image/InterpolationAdaptorInterface.hpp>
#include <core/utilities.hpp>


namespace imaging
{
  namespace polynomial_interpolation_adaptor_impl
  {
    template <class image_t>
    typename image_t::data_t interpolate_1d(const image_t & image_reference, const ublas::fixed_vector<size_t, 1> & base_coordinates, const ublas::fixed_vector<float_t, 1> & truncated_coordinates)
    {
      return typename image_t::data_t(0.0); // TODO: implement function
    }
    
    template <class image_t>
    typename image_t::data_t interpolate_2d(const image_t & image_reference, const ublas::fixed_vector<size_t, 2> & base_coordinates, const ublas::fixed_vector<float_t, 2> & truncated_coordinates)
    {
      return (1.0 - truncated_coordinates(0)) * (1.0 - truncated_coordinates(1)) * image_reference[base_coordinates] +
             truncated_coordinates(0) * (1.0 - truncated_coordinates(1)) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 0)] + 
             (1.0 - truncated_coordinates(0)) * truncated_coordinates(1) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(0, 1)] +
             truncated_coordinates(0) * truncated_coordinates(1) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 1)];   
    }
    
    
    template <class image_t>
    typename image_t::data_t interpolate_3d(const image_t & image_reference, const ublas::fixed_vector<size_t, 3> & base_coordinates, const ublas::fixed_vector<float_t, 3> & truncated_coordinates)
    {
      float_t temp(0.0);
      for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 2; ++j)
          for(size_t k = 0; k < 2; ++k)
          { 
            temp += (delta(1, i) * truncated_coordinates(0) + delta(0, i) * (1.0 - truncated_coordinates(0))) *
                    (delta(1, j) * truncated_coordinates(1) + delta(0, j) * (1.0 - truncated_coordinates(1))) *
                    (delta(1, k) * truncated_coordinates(2) + delta(0, k) * (1.0 - truncated_coordinates(2))) * 
                    image_reference[base_coordinates + ublas::fixed_vector<size_t, 3>(i, j, k)];  
          }
          
      return temp;
    }
    
    template <class image_t, class result_t>
    void interpolate_gradient_1d(const image_t & image_reference, const ublas::fixed_vector<size_t, 1> & base_coordinates, const ublas::fixed_vector<float_t, 1> & truncated_coordinates, result_t & value)
    {}
    
    template <class image_t, class result_t>
    void interpolate_gradient_2d(const image_t & image_reference, const ublas::fixed_vector<size_t, 2> & base_coordinates, const ublas::fixed_vector<float_t, 2> & truncated_coordinates, result_t & value)
    {  
      value(0) = - (1.0 - truncated_coordinates(1)) * image_reference[base_coordinates]
                 + (1.0 - truncated_coordinates(1)) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 0)]
                 - truncated_coordinates(1) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(0, 1)]
                 + truncated_coordinates(1) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 1)];
      
      value(1) = - (1.0 - truncated_coordinates(0)) * image_reference[base_coordinates]
                 - truncated_coordinates(0) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 0)]
                 + (1.0 - truncated_coordinates(0)) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(0, 1)]
                 + truncated_coordinates(0) * image_reference[base_coordinates + ublas::fixed_vector<size_t, 2>(1, 1)];
    }
    
    template <class image_t, class result_t>
    void interpolate_gradient_3d(const image_t & image_reference, const ublas::fixed_vector<size_t, 3> & base_coordinates, const ublas::fixed_vector<float_t, 3> & truncated_coordinates, result_t & value)
    {
      value.assign(0);
      
      for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 2; ++j)
          for(size_t k = 0; k < 2; ++k)
          { 
            value(0) += (delta(1, i)-delta(0, i)) *
            (delta(1, j)*truncated_coordinates(1)+delta(0, j)*(1-truncated_coordinates(1))) *
            (delta(1, k)*truncated_coordinates(2)+delta(0, k)*(1-truncated_coordinates(2))) * 
            image_reference[base_coordinates + ublas::fixed_vector<size_t, 3>(i, j, k)]; 
            
            value(1) += (delta(1, i)*truncated_coordinates(0)+delta(0, i)*(1-truncated_coordinates(0))) *
            (delta(1, j)-delta(0, j)) *
            (delta(1, k)*truncated_coordinates(2)+delta(0, k)*(1-truncated_coordinates(2))) * 
            image_reference[base_coordinates + ublas::fixed_vector<size_t, 3>(i, j, k)];  
            
            value(2) += (delta(1, i)*truncated_coordinates(0)+delta(0, i)*(1-truncated_coordinates(0))) *
            (delta(1, j)*truncated_coordinates(1)+delta(0, j)*(1-truncated_coordinates(1))) *
            (delta(1, k)-delta(0, k)) * 
            image_reference[base_coordinates + ublas::fixed_vector<size_t, 3>(i, j, k)];  
          }
    } 
  }
  
  
  /** \ingroup image
      \brief Polynomial interpolation of image values at floating point coordinates.
      
      This class provides a pixel accessor [] for real valued pixel coordinates via multilinear interpolation of the surrounding pixels.
      For planar images this mean bilinear interpolation and for volume data trilinear interpolation.
  */
  template <class image_t>
  class PolynomialInterpolationAdaptor : public InterpolationAdaptorInterface<image_t>
  {
  public:
    /** Constructs a PolynomialInterpolationAdaptor object from \em image_reference. The class \em image_t must implement ImageInterface. */
    PolynomialInterpolationAdaptor(const image_t & image_reference) : InterpolationAdaptorInterface<image_t>(image_reference)
    {}
    
    typename InterpolationAdaptorInterface<image_t>::data_t operator[](const ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> & index)
    {
      ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> local_index;
      if( ! interpolation_adaptor_impl::offset_position(InterpolationAdaptorInterface<image_t>::size(), index, local_index))
      {
        return typename InterpolationAdaptorInterface<image_t>::data_t(0);
      }
      
      ublas::fixed_vector<size_t, InterpolationAdaptorInterface<image_t>::dimension> base_coordinates;
      ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> truncated_coordinates;
      
      interpolation_adaptor_impl::compute_local_coordinates(local_index, InterpolationAdaptorInterface<image_t>::_image_reference.size(), base_coordinates, truncated_coordinates);
          
      switch(InterpolationAdaptorInterface<image_t>::dimension)
      {
        case 1:
          return polynomial_interpolation_adaptor_impl::interpolate_1d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates);
        case 2:
          return polynomial_interpolation_adaptor_impl::interpolate_2d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates);
        case 3:
          return polynomial_interpolation_adaptor_impl::interpolate_3d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates);
          
        default:
          throw Exception("Exception: PolynomialInterpolationAdaptor not implemented for image dimension " + boost::lexical_cast<std::string>(InterpolationAdaptorInterface<image_t>::dimension) + ".");
      }
    }
    
    /** Computes the gradient of the interpolation at the position \em index of the surrounding pixel values and stores it in \em value. */
    void gradient (const ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> & index, ublas::fixed_vector<typename InterpolationAdaptorInterface<image_t>::data_t, InterpolationAdaptorInterface<image_t>::dimension> & value)
    {
      typedef ublas::fixed_vector<typename InterpolationAdaptorInterface<image_t>::data_t, InterpolationAdaptorInterface<image_t>::dimension> return_data_t;
      
      ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> local_index;
      if( ! interpolation_adaptor_impl::offset_position(InterpolationAdaptorInterface<image_t>::size(), index, local_index))
      {
        value = return_data_t(0);
        return;
      }
      
      ublas::fixed_vector<size_t, InterpolationAdaptorInterface<image_t>::dimension> base_coordinates;
      ublas::fixed_vector<float_t, InterpolationAdaptorInterface<image_t>::dimension> truncated_coordinates;
      
      interpolation_adaptor_impl::compute_local_coordinates(local_index, InterpolationAdaptorInterface<image_t>::_image_reference.size(), base_coordinates, truncated_coordinates);
      
      switch(InterpolationAdaptorInterface<image_t>::dimension)
      {
        case 1:
          polynomial_interpolation_adaptor_impl::interpolate_gradient_1d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates, value);
          break;
        case 2:
          polynomial_interpolation_adaptor_impl::interpolate_gradient_2d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates, value);
          break;
        case 3:
          polynomial_interpolation_adaptor_impl::interpolate_gradient_3d(InterpolationAdaptorInterface<image_t>::_image_reference, base_coordinates, truncated_coordinates, value);
          break;
          
        default:
          throw Exception("Exception: LinearInterpolationAdaptor not implemented for image dimension " + boost::lexical_cast<std::string>(InterpolationAdaptorInterface<image_t>::dimension) + ".");
      }
    }
  };
}

#endif

