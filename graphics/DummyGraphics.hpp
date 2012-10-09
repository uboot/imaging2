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

#ifndef GRAPHICS_DUMMYGRAPHICS_H
#define GRAPHICS_DUMMYGRAPHICS_H

#include <graphics/GraphicsInterface.hpp>
#include <image/Color.hpp>
#include <image/Image.hpp>

namespace imaging
{
  /** \ingroup graphics 
      \brief Graphics output that ignores any input.
  */
  class DummyGraphics : public GraphicsInterface
  {
    static const StreamStatus stream_status;
    
  public:
    /** Standard graphics output. Currently only output to \em out is possible. */
    static DummyGraphics out;

    /** Image formats to be used in write_image(). */
    enum image_formats { 
      IMAGE_FORMAT_EPS /** EPS file with suffix \em eps. */, 
      IMAGE_FORMAT_SVG /** Compressed SVG file with suffix \em svgz. */,
      IMAGE_FORMAT_PDF /** PDF file with suffix \em pdf. */,
      IMAGE_FORMAT_DETERMINE /** Determine image format from file name suffix. */,
      /** \cond */
      IMAGE_FORMAT_CAN_NOT_DETERMINE
      /** \endcond */
    };
    
    static const size_t STD_WINDOW_SIZE = 400;

    DummyGraphics() {}
    
    virtual ~DummyGraphics() {}

    /** Initializes graphics. The arguments \em lower_left and \em upper_right specify the sections of the x-axis and the y-axis which are displayed, \em window_width and \em window_height the size of the output window. */
    virtual void init(int argc, char** argv,
                               const ublas::fixed_vector<float_t, 2> & lower_left,
                               const ublas::fixed_vector<float_t, 2> & upper_right,
                               const size_t window_width = STD_WINDOW_SIZE,
                               const size_t window_height = STD_WINDOW_SIZE) {}
     
    virtual void clear() {}
                              
    virtual void set_coordinates(const ublas::fixed_vector<float_t, 2> & lower_left, const ublas::fixed_vector<float_t, 2> & upper_right) {}

    virtual void circle(const ublas::fixed_vector<float_t, 2> & center, float_t radius) {}

    virtual void polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) {}
    
    virtual void fill_polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) {}
    
    virtual void polyline(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) {}
    
    virtual void vertex(const ublas::fixed_vector<float_t, 2> & vertex) {}
    
    virtual void image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> x_interval, const ublas::fixed_vector<float_t, 2> y_interval) {}
    
    virtual void spline_curve(const Bspline< ublas::fixed_vector<float_t, 2> > & spline_curve) {}

    virtual GraphicsInterface & operator<<(const CommandInterface & command) { return *this; }
    
    virtual GraphicsInterface & operator<<(const StreamStatus & status) { return *this; }

    virtual const StreamStatus & get_stream_status() { return stream_status; }

    /** Write the content of the graphics output to an image file. The file format of the output image is determined by \em image_format. The size of the image will be the current size of the coordinate axes. */
    virtual void  write_image(const std::string & file_name,
                      const int image_format = IMAGE_FORMAT_DETERMINE) {}

    /** Write the content of the graphics output to an image file. The file format of the output image is determined by \em image_format. The size of the image is determined by \em width and \em height. */
    virtual void  write_image(const std::string & file_name, size_t width, size_t height,
                      const int image_format = IMAGE_FORMAT_DETERMINE) {}
  };
}

#endif
