#include <graphics/OpenGlViewer.hpp>

#include <graphics/Objects.hpp>

namespace imaging
{
  namespace open_gl_viewer_impl
  {
    void corners_to_intervals(const ublas::fixed_vector<float_t, 2> & lower_left, 
                              const ublas::fixed_vector<float_t, 2> & upper_right, 
                              ublas::fixed_vector<float_t, 2> & x_axis_interval, 
                              ublas::fixed_vector<float_t, 2> & y_axis_interval)
    {
      x_axis_interval(0) = lower_left(0);
      x_axis_interval(1) = upper_right(0);
    
      y_axis_interval(0) = lower_left(1);
      y_axis_interval(1) = upper_right(1);
    }
  }
  
  const float_t OpenGlViewer::PI = 3.1415926;
  OpenGlViewer OpenGlViewer::out;
  int OpenGlViewer::_argc = 0;
  char** OpenGlViewer::_argv = 0;
  size_t OpenGlViewer::_window_width = 300;
  size_t OpenGlViewer::_window_height = 300;
  bool OpenGlViewer::_glut_active = false;
  
  void OpenGlViewer::lock_objects()
  {
    if (pthread_mutex_lock (&_object_mutex) != 0)
      throw Exception("Exception: Could not lock pthread mutex.");
  }

  void OpenGlViewer::unlock_objects()
  {
    if (pthread_mutex_unlock (&_object_mutex) != 0)
      throw Exception("Exception: Could not unlock pthread mutex.");
  }

  OpenGlViewer::OpenGlViewer() :
    _redisplay(0), 
    _current_stream_status(),
    _current_object_index(0), 
    _x_axis_interval(0.0, 1.0),
    _y_axis_interval(0.0, 1.0),
    _exit_mode(false)
  {
    if( pthread_mutex_init(&_object_mutex, NULL) != 0 )
      throw Exception("Exception: Could not create mutex.");
    
    *this << StreamStatus();
  }

  OpenGlViewer::~OpenGlViewer()
  {
    pthread_mutex_destroy( &_object_mutex );
  }

  void OpenGlViewer::init(int argc, char** argv,
                          const ublas::fixed_vector<float_t, 2> & lower_left,
                          const ublas::fixed_vector<float_t, 2> & upper_right,
                          size_t window_width,
                          size_t window_height)
  {
    lock_objects();
    open_gl_viewer_impl::corners_to_intervals(lower_left, upper_right, _x_axis_interval, _y_axis_interval);
    _argc = argc;
    _argv = argv;
    _window_width = window_width;
    _window_height = window_height;
    unlock_objects();
  }
  
  void * OpenGlViewer::display(void *)
  {
    glutInit(&_argc, _argv);

    glutInitDisplayMode(GLUT_SINGLE | GLUT_DEPTH);
    glutInitWindowSize(int(_window_width), int(_window_height));
    glutInitWindowPosition(100,100);
    glutCreateWindow(_argv[0]);
    _glut_active = true;

    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_POLYGON_OFFSET_LINE);
    glDepthFunc(GL_LESS);
    //     glEnable(GL_LINE_SMOOTH);
    //     glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    
    glutDisplayFunc(&OpenGlViewer::static_display_function);
    glutKeyboardFunc(&OpenGlViewer::static_keyfunc);
    glutReshapeFunc(&OpenGlViewer::static_reshape);
    glutIdleFunc(&OpenGlViewer::static_idle_function);
    
    glutMainLoop();
    
    _glut_active = false;
    
    return 0;
  }
  
  void OpenGlViewer::idle_function()
  {
    lock_objects();
    if(_redisplay > 0)
    {
      glutPostRedisplay();
      _redisplay--;
    }
    unlock_objects();
    
    usleep(100000);
  }
  
  void OpenGlViewer::reshape(int w, int h)
  {
    lock_objects();
    set_viewport(w, h);
    unlock_objects();
  }
  
  void OpenGlViewer::display_function()
  {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    lock_objects();
    
    std::set< boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> >::iterator new_object_iter;
    for(new_object_iter = _new_objects.begin(); new_object_iter != _new_objects.end(); ++new_object_iter)
      (*new_object_iter)->initialize();
    _new_objects.clear();
    
    std::size_t clear_index = 0;
    std::map<size_t, boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> >::iterator object_iter;
    for(object_iter = _objects.begin(); object_iter != _objects.end(); ++object_iter)
      if(dynamic_cast<open_gl_viewer_impl::ClearPipeline * >(object_iter->second.get()))
        clear_index = object_iter->first;
      else
        object_iter->second->execute(open_gl_viewer_impl::SCREEN_MODE);
    
    _objects.erase(_objects.begin(), _objects.find(clear_index));
      
    unlock_objects();
    
    glFlush();
  }

  void OpenGlViewer::keyfunc(unsigned char key, int x, int y)
  {
    switch (int(key))
    {
    case RETURN_LINUX:
    case ESC:
      _exit_mode = true;
      exit(0);
      break;
    default:;
    }
  }
  
  void OpenGlViewer::clear()
  {
    lock_objects();
    _redisplay = 0;
    _current_stream_status = StreamStatus();
    // _current_object_index = 0;
    _x_axis_interval = ublas::fixed_vector<float_t, 2>(0.0, 1.0);
    _y_axis_interval = ublas::fixed_vector<float_t, 2>(0.0, 1.0);
    
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::ClearPipeline()));
    
    unlock_objects();
    
    *this << StreamStatus();
    
    lock_objects();
    request_redisplay();
    unlock_objects();
  }
  
  void OpenGlViewer::request_redisplay()
  {
    _redisplay++;
  }
  
  void OpenGlViewer::add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> object)
  {
    if(_exit_mode)
      return;
      
    _objects[_current_object_index++] = object;
    _new_objects.insert(object);
  }
  
  void OpenGlViewer::add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> object)
  {
    if(_exit_mode)
      return;
      
    _objects[_current_object_index] = object;
    _new_objects.insert(object);
    _group_object_map.insert(std::pair<size_t, size_t>(_current_stream_status.group_id(), _current_object_index));
    _current_object_index++;
  }
                               
  void OpenGlViewer::set_coordinates(const ublas::fixed_vector<float_t, 2> & lower_left, const ublas::fixed_vector<float_t, 2> & upper_right)
  {
    lock_objects();
    
    open_gl_viewer_impl::corners_to_intervals(lower_left, upper_right, _x_axis_interval, _y_axis_interval);   
    
    if(_glut_active)
    { 
      int width = glutGet(GLUT_WINDOW_WIDTH);
      int height = glutGet(GLUT_WINDOW_HEIGHT);
      glutReshapeWindow(width, height);
    }
      
    unlock_objects();
  }

  void OpenGlViewer::circle(const ublas::fixed_vector<float_t, 2> & center, float_t radius)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > vertices(NUM_CIRCLE_POINTS);

    for(size_t i = 0; i < NUM_CIRCLE_POINTS; i++)
      vertices[i] = center + radius * ublas::fixed_vector<float_t, 2>(cos(i * 2 * PI / NUM_CIRCLE_POINTS),
                    sin(i * 2 * PI / NUM_CIRCLE_POINTS));
    polygon(vertices);
  }
  
  void OpenGlViewer::polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::Polygon(vertices)));
    unlock_objects();
  }
  
  void OpenGlViewer::fill_polygon(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::FillPolygon(vertices)));
    unlock_objects();
  }
  
  void OpenGlViewer::polyline(const std::vector< ublas::fixed_vector<float_t, 2> > & vertices)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::Polyline(vertices)));
    unlock_objects();
  }
  
  void OpenGlViewer::vertex(const ublas::fixed_vector<float_t, 2> & vertex)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::Vertex(vertex)));
    unlock_objects();
  }
  
  void OpenGlViewer::image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> x_interval, const ublas::fixed_vector<float_t, 2> y_interval)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::Image(image, x_interval, y_interval)));
    unlock_objects();
  }
  
  void OpenGlViewer::spline_curve(const Bspline< ublas::fixed_vector<float_t, 2> > & spline_curve)
  {
    lock_objects();
    add_to_group(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::SplineCurve(spline_curve)));
    unlock_objects();
  }

  GraphicsInterface & OpenGlViewer::operator<<(const StreamStatus & status)
  {
    lock_objects();
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::SetColor(status.color())));
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new
      open_gl_viewer_impl::OffsetZLayer(status.offset() - _current_stream_status.offset())));
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new
      open_gl_viewer_impl::Translate(status.translation() - _current_stream_status.translation())));
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::LineWidth(status.line_width())));

    _current_stream_status = status;
    unlock_objects();
    return *this;
  }
  
  GraphicsInterface & OpenGlViewer::operator<<(const CommandInterface & command)
  {
    lock_objects();
    
    if(dynamic_cast<const FlushCommand *>(&command))
      request_redisplay();
      
    if(dynamic_cast<const ResetCommand *>(&command))
    {
      while(_group_object_map.count(_current_stream_status.group_id()))
      {
        std::multimap<size_t, size_t>::iterator iter = _group_object_map.find(_current_stream_status.group_id());
        if(open_gl_viewer_impl::Image * image = dynamic_cast<open_gl_viewer_impl::Image *>(_objects[iter->second].get()))
          add_to_stream(boost::shared_ptr<open_gl_viewer_impl::DeleteTexture>(new open_gl_viewer_impl::DeleteTexture(image->texture_name())));
        _objects.erase(iter->second);
        _group_object_map.erase(iter);
      }
      
      request_redisplay();
    }
      
    if(const GraphicsInterface::set_line_width * cmd = dynamic_cast<const GraphicsInterface::set_line_width *>(&command))
    {
      add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::LineWidth(cmd->line_width())));
      _current_stream_status.line_width(cmd->line_width());
    }
      
    if(const GraphicsInterface::set_color * cmd = dynamic_cast<const GraphicsInterface::set_color *>(&command))
    {
      add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::SetColor(cmd->color())));
      _current_stream_status.color(cmd->color());
    }
      
    if(const GraphicsInterface::offset_z_layer * cmd = dynamic_cast<const GraphicsInterface::offset_z_layer *>(&command))
    {
      add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::OffsetZLayer(cmd->offset())));
      _current_stream_status.offset(_current_stream_status.offset() + cmd->offset());
    }
      
    if(const GraphicsInterface::translate * cmd = dynamic_cast<const GraphicsInterface::translate *>(&command))
    {
      add_to_stream(boost::shared_ptr<open_gl_viewer_impl::ObjectInterface>(new open_gl_viewer_impl::Translate(cmd->translation())));
      _current_stream_status.translation(_current_stream_status.translation() + cmd->translation());
    }
      
    if(const GraphicsInterface::set_group * cmd = dynamic_cast<const GraphicsInterface::set_group *>(&command))
    {
      _current_stream_status.group_id(cmd->group_id());
      request_redisplay();
    }
    
    unlock_objects();
    
    return *this;
  }

  void OpenGlViewer::write_image(const std::string & file_name, size_t width, size_t height, const int image_format)
  {
    lock_objects();
    add_to_stream(boost::shared_ptr<open_gl_viewer_impl::WriteImage>(new open_gl_viewer_impl::WriteImage(file_name, width, height, _x_axis_interval, _y_axis_interval, image_format, out)));
    unlock_objects();
  }
  
  void OpenGlViewer::set_viewport(int width, int height)
  {
    set_viewport(width, height, _x_axis_interval, _y_axis_interval);
  }
    
  void OpenGlViewer::set_viewport(int width, int height, const ublas::fixed_vector<float_t, 2> & x_axis_interval, const ublas::fixed_vector<float_t, 2> & y_axis_interval)
  {
    GLfloat wh_ratio_window = 1.0;
    if ( ( height > 0) && ( width > 0) )
      wh_ratio_window = GLfloat(width) / GLfloat(height);

    GLfloat wh_ratio_coordinates = 1.0;
    if ( ( x_axis_interval(1) - x_axis_interval(0) > 0 ) &&
         ( y_axis_interval(1) - y_axis_interval(0) > 0 ) )
      wh_ratio_coordinates = GLfloat( x_axis_interval(1) - x_axis_interval(0) ) /
                             GLfloat( y_axis_interval(1) - y_axis_interval(0) );


    glMatrixMode(GL_PROJECTION);
    glViewport(0, 0, width, height);
    glLoadIdentity();
    if (wh_ratio_window > wh_ratio_coordinates)
      // window wider than coordinate system
      gluOrtho2D(wh_ratio_window / wh_ratio_coordinates * x_axis_interval(0),
                 wh_ratio_window / wh_ratio_coordinates * x_axis_interval(1),
                 y_axis_interval(0), y_axis_interval(1));
    else
      // window higher than coordinate system
      gluOrtho2D(x_axis_interval(0), x_axis_interval(1),
                 wh_ratio_coordinates / wh_ratio_window * y_axis_interval(0),
                 wh_ratio_coordinates / wh_ratio_window * y_axis_interval(1));

    glMatrixMode(GL_MODELVIEW);
  }
}

