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

#ifndef GRAPHICS_OPENGLVIEWER_H
#define GRAPHICS_OPENGLVIEWER_H

#include <graphics/ObjectInterface.hpp>
#include <graphics/DummyGraphics.hpp>

#include <map>
#include <set>
#include <boost/shared_ptr.hpp>
#include <pthread.h>


namespace imaging
{
  namespace open_gl_viewer_impl
  {
    class WriteImage;
    class Image;
  }

  /** \ingroup graphics
      \brief Graphics output that renders graphics using GLUT.
  */
  class OpenGlViewer : public DummyGraphics
  {
    friend class open_gl_viewer_impl::WriteImage;
    
    static const float_t PI;
    static const int NUM_CIRCLE_POINTS = 64;

    pthread_mutex_t _object_mutex;

    size_t _redisplay;
    StreamStatus _current_stream_status;
    
    std::map<size_t, boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> > _objects;
    std::set<boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> > _new_objects; // objects which have not been initialized
    std::multimap<size_t, size_t> _group_object_map;
    size_t _current_object_index;
    bool _exit_mode;
    
    static int _argc;
    static char** _argv;
    static size_t _window_width;
    static size_t _window_height;
    static bool _glut_active;
    
    ublas::fixed_vector<float_t, 2> _x_axis_interval;
    ublas::fixed_vector<float_t, 2> _y_axis_interval;
    
    enum keycodes { ESC = 27, RETURN_LINUX = 13 };
    
    /* This functions assume that the objects are locked. BEGIN */
    /* GLUT thread only BEGIN */
    static void static_idle_function() { out.idle_function(); }
    static void static_reshape(int w, int h) { out.reshape(w, h); }
    static void static_display_function() { out.display_function(); }
    static void static_keyfunc(unsigned char key, int x, int y) { out.keyfunc(key, x, y); }
    
    void idle_function();
    void reshape(int w, int h);
    void display_function();
    void keyfunc(unsigned char key, int x, int y);
    /* GLUT thread only END */
    
    void request_redisplay();
    void set_viewport(int width, int height);
    void set_viewport(int width, int height, const ublas::fixed_vector<float_t, 2> & x_axis_interval, const ublas::fixed_vector<float_t, 2> & y_axis_interval);
    
    void add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> object);
    void add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> object);
    /* This functions assume that the objects are locked. END */
    
    void lock_objects();
    void unlock_objects();

  public:
    static OpenGlViewer out;
    
    static void * display(void *);

    OpenGlViewer();

    virtual ~OpenGlViewer();

    virtual void init(int argc, char** argv,
                               const ublas::fixed_vector<float_t, 2> & lower_left,
                               const ublas::fixed_vector<float_t, 2> & upper_right,
                               size_t window_width = STD_WINDOW_SIZE,
                               size_t window_height = STD_WINDOW_SIZE);
                               
    virtual void clear();
                               
    virtual void set_coordinates(const ublas::fixed_vector<float_t, 2> & lower_left, const ublas::fixed_vector<float_t, 2> & upper_right);

    virtual void circle(const ublas::fixed_vector<float_t, 2> & center, float_t radius);
    
    virtual void polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices);
    
    virtual void fill_polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices);
    
    virtual void polyline(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices);
    
    virtual void vertex(const ublas::fixed_vector<float_t, 2> & vertex);
    
    virtual void image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> x_interval, const ublas::fixed_vector<float_t, 2> y_interval);
    
    virtual void spline_curve(const Bspline< ublas::fixed_vector<float_t, 2> > & spline_curve);

    virtual GraphicsInterface & operator<<(const StreamStatus & status);
    
    virtual GraphicsInterface & operator<<(const CommandInterface & command);

    virtual const StreamStatus & get_stream_status() { return _current_stream_status; }

    virtual void write_image(const std::string & file_name, const int image_format = IMAGE_FORMAT_DETERMINE)
    {
      write_image(file_name,
                  size_t(_x_axis_interval(1) - _x_axis_interval(0) + 0.5),
                  size_t(_y_axis_interval(1) - _y_axis_interval(0) + 0.5),
                  image_format);
    }
    
    virtual void write_image(const std::string & file_name, size_t width, size_t height,
                     const int image_format = IMAGE_FORMAT_DETERMINE);
  };
}

#endif

