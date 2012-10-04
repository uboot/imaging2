// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_SCALARIMAGE_H
#define IMAGE_SCALARIMAGE_H

#include <image/ImageInterface.hpp>

namespace imaging
{

  /** \ingroup image
      \brief %Image with constant pixel values.
      
      This class models images with constant pixel values, i.e. every pixel has the same value. This value is set during construction. The size (in memory) of ScalarImage objects is independent of their size (as an image). 
  */
  template <std::size_t N, class DATA_t>
  class ScalarImage : public ImageInterface<N, DATA_t>
  {
  public:
    static const std::size_t dimension = N;
    typedef DATA_t data_t;

  private:
    ublas::fixed_vector<size_t, dimension> _size;
    data_t _value;

  public:
    /** Constructs a scalar image of size \em size. Each pixel of this image will evaluate to \em value. */
    ScalarImage(const ublas::fixed_vector<size_t, dimension> & size, data_t value) :
      _size(size), _value(value) {}
      
    /** Copy constructor. The complexity of this is function is constant (i.e. independent of the image size). */
    ScalarImage(const ScalarImage<dimension, data_t> & source) :
      _size(source._size), _value(source._value) {}
    
    const DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> &) const { return _value; }
    
    const ublas::fixed_vector<size_t, dimension> & size() const { return _size; }
  };


}

#endif
