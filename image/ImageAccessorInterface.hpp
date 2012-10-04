// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef IMAGE_IMAGEACCESSORINTERFACE_H
#define IMAGE_IMAGEACCESSORINTERFACE_H

#include <image/ImageInterface.hpp>

namespace imaging
{

  /** \ingroup image
      \brief Abstract base class of all image accessor classes.
      
      This class defines an interface for <em>image accessor classes</em>. These are classes which provide layered access to objects which implement ImageInterface. E.g. such a class might convert values of a colour image to grayscale and back again in a transparent fashion when its pixel accessor functions (i.e. %operator[]) are called. Similar tasks include rescaling image values or selecting channels in vector valued data.
      
      Note that objects implementing this interface must also implement ImageInterface.
      
      If you want to construct an accessor from constant images, it is necessary that the ImageAccessorInterface template is instantiated for the \em const image data type:
  \code
  const GrayImage2d image;
  
  CastAccessor<const GrayImage2d, float_t> image_accessor;  // works fine
  CastAccessor<GrayImage2d, float_t> image_accessor; // will not compile!
  \endcode
      
      \sa CastAccessor, ScalingAccessor, VectorComponentAccessor
  */
  template <class image_t, class DATA_t>
  class ImageAccessorInterface : public ImageInterface<image_t::dimension, DATA_t>
  {
  protected:
    /** \brief Reference to the image object interfaced by this object. */
    image_t & _image_reference;

  public:
    static const size_t dimension = image_t::dimension;
    typedef DATA_t data_t;

    /** Initialize the image accessor with an reference to \em image. By \em image_t we mean any type, which implements ImageInterface, i.e. \em image may itself be an image accessor. This reference will be stored and accessed when the image accessor is utilized. The dimensions of the accessor are the same as those of \em image.
    */
    ImageAccessorInterface(image_t & image) : _image_reference(image) {}

    const ublas::fixed_vector<size_t, dimension> & size() const { return _image_reference.size(); }
    
    /** Returns true if the image is empty, i.e. if at least on of its dimensions is zero. */
    bool empty() const
    {
      for(size_t i = 0; i < dimension; ++i)
        if (size()(i) == 0)
          return true;
      
      return false;
    }
  };

}

#endif
