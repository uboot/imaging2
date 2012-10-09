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

#ifndef GRAPHICS_GRAPHICSINTERFACE_H
#define GRAPHICS_GRAPHICSINTERFACE_H

#include <image/Color.hpp>
#include <image/Image.hpp>
#include <spline/Bspline.hpp>

namespace imaging
{
  class OpenGlViewer;
  class DummyGraphics;
  
  /** \ingroup graphics 
      \brief Abstract class interface for graphics output.
  */
  class GraphicsInterface
  {
  public:
  
    /** \brief Graphic stream commands. */
    class CommandInterface
    {
    public:
      virtual ~CommandInterface() {}
    };
    
    class FlushCommand : public CommandInterface
    {};
    
    class ResetCommand : public CommandInterface
    {};

    /** \brief Graphics command to set the current group. */
    class set_group : public CommandInterface
    {
      size_t _group_id;
      
    public:
      /** Passing this command to a graphics stream sets the current group to \em group_id. */
      set_group(size_t group_id) : _group_id(group_id) {}
      
      /** Returns the group. */
      size_t group_id() const { return _group_id; }
      
      /** Sets the group. */
      void group_id(size_t id) { _group_id = id; }
    };


    /** \brief Graphics command to offset the current z-value. */
    class offset_z_layer : public CommandInterface
    {
      int _offset;
      
    public:
      /** \brief Passing this command to a graphics stream offsets the current z-value by \em offset. Negative values move the active layer to the back and positive values to the front.
      */
      offset_z_layer(int offset) : _offset(offset) {}
      
      /**  Returns the offset. */
      int offset() const { return _offset; }
      
      /**  Sets the offset. */
      void offset(size_t offset) { _offset = offset; }
    };

    /** \brief Graphics command to translate the current coordinate system. */
    class translate : public CommandInterface
    {
      friend class GraphicsInterface;
      
      ublas::fixed_vector<float_t, 2> _translation;
    public:
      /** \brief Passing this command to a graphics stream moves the current coordinate system by \em offset. */
      translate(const ublas::fixed_vector<float_t, 2> & translation) : _translation(translation) {}
      /** Returns the translation. */
      const ublas::fixed_vector<float_t, 2> & translation() const { return _translation; }
      
      /** Sets the translation. */
      void translation(const ublas::fixed_vector<float_t, 2> & translation) { _translation = translation; }
    };
    
    /** \brief Graphics command to set the current color. */
    class set_color : public CommandInterface
    {
      friend class GraphicsInterface;
      
      Color _color;
    public:
      /** Passing this command to a graphics stream sets the current color to \em color. */
      set_color(const Color & color) : _color(color) {}
      
      /** Returns the color. */
      const Color & color() const { return _color; }
      
      /** Sets the color. */
      void color(const Color & color) { _color = color; }
    };

    /** \brief Graphics command to set the current line width. */
    class set_line_width : public CommandInterface
    {
      friend class GraphicsInterface;
      
      float_t _line_width;
    public:
      /** Passing this command to a graphics stream sets the current line width to \em line_width. */
      set_line_width(float_t line_width) : _line_width(line_width) {}
      
      /** Returns the line width. */
      float_t line_width() const { return _line_width; }
      
      /** Sets the line width. */
      void line_width(float_t width) { _line_width = width; }
    };
    
       /** \brief Class which stores the current status of the graphics output. */
    class StreamStatus : public set_group, public offset_z_layer, public set_color, public set_line_width, public translate
    {

    public:
      StreamStatus() : set_group(0), offset_z_layer(0), set_color(Color::BLACK), set_line_width(1.0), translate(ublas::fixed_vector<float_t, 2>(0.0, 0.0)) {}
      
      /** Constructs and initializes a stream status. */
      StreamStatus(size_t group_id, int offset, const Color & color, float_t line_width, const ublas::fixed_vector<float_t, 2> & translation) : set_group(group_id), offset_z_layer(offset), set_color(color), set_line_width(line_width), translate(translation) {}
    };
    
  public:
    
    /** Graphics command to flush the current group. */
    static const FlushCommand flush;
    
    /** Graphics command to reset the current group. */
    static const ResetCommand reset_group;
    
    virtual ~GraphicsInterface() {}
     
    /** Clears the graphics output. All groups are reset and all objects removed from the window. */                          
    virtual void clear() = 0;
                              
    /** Initializes the coordinate system of the graphics output. The arguments \em lower_left and \em upper_right specify the sections of the x-axis and the y-axis which are displayed. */
    virtual void set_coordinates(const ublas::fixed_vector<float_t, 2> & lower_left, const ublas::fixed_vector<float_t, 2> & upper_right) = 0;

    /** Draws the line segment determined by \em v_0 and \em v_1. */
    void line_segment(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1);
    
    /** Draws the vector \em direction at \em starting_point and adds an arrowhead. */ 
    void arrow(const ublas::fixed_vector<float_t, 2> & starting_point, const ublas::fixed_vector<float_t, 2> & direction) ;
    
    /** Draws the quadrangle determined by \em v_0, \em v_1, \em v_2 and \em v_3. */
    void quadrangle(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1, const ublas::fixed_vector<float_t, 2> & v_2, const ublas::fixed_vector<float_t, 2> & v_3);
    
    /** Draws the triangle determined by \em v_0, \em v_1 and \em v_2. */
    void triangle(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1, const ublas::fixed_vector<float_t, 2> & v_2);
    
    /** Draws the circle determined by \em center and \em radius.
        \sa operator<<(gr & out, const Circle & circle)
    */
    virtual void circle(const ublas::fixed_vector<float_t, 2> & center, float_t radius) = 0;

    /** Draws the polygon determined by \em vertices. The first and the last point in \em vertices will automatically be connected, i.e. there is no need for them to be the same.
    \sa polyline()*/
    virtual void polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) = 0;

    /** Fills the polygon determined by \em vertices with the current color. */
    virtual void fill_polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) = 0;
    
    /** Draws the polyline determined by \em vertices. The first and the last point in \em vertices will not be connected.
    \sa polygon */
    virtual void polyline(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices) = 0;
   
    /** Draws \em vertex as point on the drawing plane.
        \sa operator<<(gr & out, const ublas::fixed_vector<float_t, 2> &)
    */
    virtual void vertex(const ublas::fixed_vector<float_t, 2> & vertex) = 0;
    
    /** Rescales \em image such that it fits into the rectangle \em x_interval x \em y_interval and draws it. 
        \sa operator<<(gr & out, const ColorImage2d &)
    */
    virtual void image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> x_interval, const ublas::fixed_vector<float_t, 2> y_interval) = 0;
    
    /** Draws the spline curve of \em order with \em knots and \em coefficient. 
        \sa operator<<(gr & out, const Bspline< ublas::fixed_vector<float_t, 2> > &), operator<<(gr & out, const PeriodicBspline< ublas::fixed_vector<float_t, 2> > &), operator<<(gr & out, const BsplineShape & spline_shape)
    */
    virtual void spline_curve(const Bspline< ublas::fixed_vector<float_t, 2> > & spline_curve) = 0;

    /** Pass the stream command \em command to the stream. */
    virtual GraphicsInterface & operator<<(const CommandInterface & command) = 0;
    
    /** Pass the stream status \em status to the stream. */
    virtual GraphicsInterface & operator<<(const StreamStatus & status) = 0;

    /** Get the current stream status. Use this function to backup the status and restore it later. The status of a stream is determined by its 
        - color,
        - group,
        - z_value,
        - translation, and
        - line width.
    */
    virtual const StreamStatus & get_stream_status() = 0;
  };
}

#endif
