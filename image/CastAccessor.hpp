#ifndef IMAGE_CASTACCESSOR_H
#define IMAGE_CASTACCESSOR_H


#include <image/ImageAccessorInterface.hpp>

namespace imaging
{

  /** \ingroup image
      \brief Provides casted pixel data type access to images.
      
      This class provides casted access to objects implementing ImageInterface. Consider e.g. a colour image that you want to use in a functions which makes sense only for scalar, i.e. grayscale images. Instead of converting the whole colour image to grayscale you can construct a CastAccessor<ColorImage2d, GrayValue> from the orginal image, which then behaves as grayscale image.
      The only prerequisite is that it is possible to explicitely cast the pixel value from the input image to \em DATA_t (and vice-versa if the user writes to the accessor). In the given example Color can be indeed casted to GrayValue and back again.
  */
  template <class image_t, class DATA_t = float_t>
  class CastAccessor : public ImageAccessorInterface<image_t, DATA_t>
  {
  public:
    static const size_t dimension = ImageAccessorInterface<image_t, DATA_t>::dimension;
    typedef DATA_t data_t;

  private:
    /** \cond */
    class data_reference
    {
      CastAccessor<image_t, DATA_t> & _accessor_reference;
      const ublas::fixed_vector<size_t, dimension> & _index_reference;

    public:
      data_reference(CastAccessor<image_t, DATA_t> & accessor, const ublas::fixed_vector<size_t, dimension> & index) : _accessor_reference(accessor), _index_reference(index) {}

      DATA_t operator=(const DATA_t & input)
      { _accessor_reference._image_reference[_index_reference] = (typename image_t::data_t)(input); return input; }

      operator DATA_t() const { return DATA_t(_accessor_reference._image_reference[_index_reference]); }
    };
    /** \endcond */

  public:
    /** Initialize the image accessor with an reference to \em image which must be an ImageInterface object. This reference will be stored and accessed when the image accessor is utilized. The dimensions of the accessor are the same as those of \em image. 
        
        Note that it must be possible to explicitely cast from \em image_t::data_t to \em DATA_t (and from \em DATA_t to \em image_t::data_t for write access) to use the accessor. */
    CastAccessor(image_t & image) : ImageAccessorInterface<image_t, DATA_t>(image) {}

    DATA_t operator[](const ublas::fixed_vector<size_t, dimension> & index) const { return DATA_t(ImageAccessorInterface<image_t, DATA_t>::_image_reference[index]); }
    data_reference operator[](const ublas::fixed_vector<size_t, dimension> & index) { return data_reference(*this, index); }
  };

}

#endif
