// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_IMAGEINTERFACE_H
#define IMAGE_IMAGEINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup image
      \brief Abstract base class of all classes modelling image data.
      
      This class defines an interface which should be implemented by classes for <em>image data</em>. By image data we mean an \em N-dimensional array. The data \em DATA_t stored in the array is arbitrary but must be the same for each entry in the array. The range of valid indices for entries of the array is fixed upon construction of the array (although derived classes, such as Image, may proved functions to resize images) and can be obtained by calling size(). I.e. range of valid entries for the <em>i</em>-th component of an image runs from 0 to size()(i) - 1. Note that the operator [] is always used the same way to access pixels, independent of the image dimensions. This is demonstrated in the following example:
      
  \code
  // float_image_2d implements ImageInterface<2, float_t>
  // graysc_image_3d implements ImageInterface<3, unsigned char>
  
  float_image_2d[ublas::fixed_vector<size_t, 2>(0, 1)] = 0.5;
  graysc_image_3d[ublas::fixed_vector<size_t, 3>(0, 1, 1)] = 128;
  \endcode
  */
  template <std::size_t N, class DATA_t>
  class ImageInterface
  {
  public:
    /** The dimension of the image. It equals 2 for planar images and 3 for voxel-data. */
    static const std::size_t dimension = N;
    
    /** The data type which is stored in the image. Common choices are \em float_t for grayscale images stored at floating point precision or \em char for 8-bit data. */
    typedef DATA_t data_t;

    /** Returns a constant reference to a pixel. */
    const DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> & index) const;
    
    /** Returns a reference to a pixel. This function might not be implemented for accessor classes to constant data, such as ConstImageAccessorInterface(). */
    DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> & index);

    /** Returns the size of the image. */
    const ublas::fixed_vector<size_t, dimension> & size() const;
  };
}

#endif
