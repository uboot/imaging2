/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace imaging
{
  /** \cond */
  template <size_t N, class DATA_t> class Image;
  
  namespace image_impl
  {
    template <size_t N>
    class image_access
    {
    public:
      template <class DATA_t>
      static const DATA_t & access(const Image<N, DATA_t> & image, const ublas::fixed_vector<size_t, N> & index);
      
      template <class DATA_t>
      static DATA_t & access(Image<N, DATA_t> & image, const ublas::fixed_vector<size_t, N> & index);
    };
    
    template <>
    class image_access<2>
    {
    public:
      template <class DATA_t>
      static const DATA_t & access(const Image<2, DATA_t> & image, const ublas::fixed_vector<size_t, 2> & index);
      
      template <class DATA_t>
      static DATA_t & access(Image<2, DATA_t> & image, const ublas::fixed_vector<size_t, 2> & index);
    };
  }
  /** \endcond */
  
  template <size_t N>
  bool increment_index(const ublas::fixed_vector<size_t, N> & dimensions, size_t fixed_component, ublas::fixed_vector<size_t, N> & index)
  {
    for(size_t i = 0; i < N; ++i)
    {
      if(i == fixed_component)
        continue;

      if(++index(i) == dimensions(i))
        if(i == N)
          return false;
        else
          index(i) = 0;
      else
        return true;
    }

    return false;
  }

  template <size_t N>
  bool increment_index(const ublas::fixed_vector<size_t, N> & dimensions, ublas::fixed_vector<size_t, N> & index)
  {
    for(size_t i = 0; i < N; ++i)
    {
      if(++index(i) == dimensions(i))
        if(i == N)
          return false;
        else
          index(i) = 0;
      else
        return true;
    }

    return false;
  }

  /** \ingroup image
      \brief Stores image data of arbitrary dimension and pixel type.
      
      Use this class template to store image data of dimension \em N with pixel data of type \em DATA_t. Objects of this class can be resized by calling set_size(). This is not generally true for implementations of ImageInterface.
  */
  template <std::size_t N, class DATA_t>
  class Image : public ImageInterface<N, DATA_t>, public boost::multi_array<DATA_t, N>
  {
  public:
    static const std::size_t dimension = N;
    typedef DATA_t data_t;

  private:
    ublas::fixed_vector<size_t, dimension> _size;

  public:
    /** Default constructor. Sets the size of every image dimension to zero. */
    Image() : boost::multi_array<DATA_t, dimension>(), _size(ublas::scalar_vector<size_t>(N, 0)) {  }

    /** Construct an image by specifying its dimensions. */
    explicit Image(const ublas::fixed_vector<size_t, dimension> & size) : boost::multi_array<DATA_t, dimension>(size), _size(size)
    { }
    
    /** Construct an image from the file \em file_name. The type of the file is determined from the suffix of \em file_name. Currently only 2-dimensional images are supported. */
    explicit Image(const std::string & file_name)
    { 
      read_image(file_name);
    }
    
    /** Construct an image from the file \em file_name. The type of the file is determined from the suffix of \em file_name. Currently only 2-dimensional images are supported. */
    explicit Image(const char * file_name)
    { 
      read_image(file_name);
    }
    
    /** Construct an image from any source which implements ImageInterface of the same dimension and whose data type can converted to the data type of this image. In other words, the explicit conversion from \em source_accessor_t::data_t to \em data_t must be defined. */
    template <class source_accessor_t>
    explicit Image(const source_accessor_t & image_accessor) : boost::multi_array<DATA_t, dimension>(image_accessor.size()), _size(image_accessor.size())
    { 
      copy_image(image_accessor, *this);
    }
    
    /** Assignement operator. */
    Image<N, DATA_t> & operator=(const Image<N, DATA_t> & source)
    { 
      resize(source.size());
      boost::multi_array<DATA_t, dimension>::operator=((const boost::multi_array<DATA_t, dimension> &)(source));
      
      return *this;
    }
    
    /** Assignement operator for any image which implements ImageInterface of the same dimension and whose data type can converted to the data type of this image. In other words, the explicit conversion from \em source_accessor_t::data_t to \em data_t must be defined. */
    template <class source_accessor_t>
    const Image<N, DATA_t> & operator=(const source_accessor_t & image_accessor)
    { 
      resize(image_accessor.size());
      copy_image(image_accessor, *this);
      
      return *this;
    }

    const DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> & index) const { return image_impl::image_access<N>::access(*this, index); }
    
    DATA_t & operator[](const ublas::fixed_vector<size_t, dimension> & index) { return image_impl::image_access<N>::access(*this, index); }

    const ublas::fixed_vector<size_t, dimension> & size() const { return _size; }
    
    /** Sets the size of the image. */
    void resize(const ublas::fixed_vector<size_t, dimension> & size)
    {
      boost::multi_array<DATA_t, dimension>::resize(size);
      _size = size;
    }
    
    /** Returns true if the image is empty, i.e. if at least on of its dimensions is zero. */
    bool empty() const
    {
      for(size_t i = 0; i < N; ++i)
        if (size()(i) == 0)
          return true;
      
      return false;
    }
    
    /** Reads the image from the file \em file_name. The type of the file is determined from the suffix of \em file_name. Currently only 2-dimensional images are supported. */
    void read_image(const std::string & file_name) { }

    /** Writes the image to the file \em file_name. The type of the file is determined from the suffix of \em file_name. Currently only 2-dimensional images are supported. */
    void write_image(const std::string & file_name) const { }
  };

  template <class source_accessor_t, class target_accessor_t>
  void copy_image(const source_accessor_t & source, target_accessor_t & target)
  {
    if(source.size() != target.size())
      throw Exception("Exception: Dimensions in copy_image() do not agree.");

    const std::size_t dimension = source_accessor_t::dimension;
    std::size_t n_pixels = 1;

    for(std::size_t i = 0; i < dimension; ++i)
      n_pixels *= source.size()(i);

    if(n_pixels == 0) return;

    ublas::fixed_vector<size_t, dimension> index;
    index.assign(0);

    do
    {
      target[index] = (typename target_accessor_t::data_t)(source[index]);
    }
    while(increment_index(source.size(), index));
  }
  
  namespace image_impl
  {
    template <size_t N>
    template <class DATA_t>
    const DATA_t & image_access<N>::access(const Image<N, DATA_t> & image, const ublas::fixed_vector<size_t, N> & index)
    {
      return image(index);
    }
    
    template <size_t N>
    template <class DATA_t>
    DATA_t & image_access<N>::access(Image<N, DATA_t> & image, const ublas::fixed_vector<size_t, N> & index)
    {
      return image(index);
    }
    
    template <class DATA_t>
    const DATA_t & image_access<2>::access(const Image<2, DATA_t> & image, const ublas::fixed_vector<size_t, 2> & index)
    {
      return ((const boost::multi_array<DATA_t, 2> &)(image))[index(0)][index(1)];
    }
    
    template <class DATA_t>
    DATA_t & image_access<2>::access(Image<2, DATA_t> & image, const ublas::fixed_vector<size_t, 2> & index)
    {
      return ((boost::multi_array<DATA_t, 2> &)(image))[index(0)][index(1)];
    }
  }
}

#endif
