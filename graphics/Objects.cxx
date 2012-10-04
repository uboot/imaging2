#include <graphics/Objects.hpp>

#include <graphics/OpenGlViewer.hpp>
#include <external/gl2ps/gl2ps.h>

namespace imaging
{
  namespace open_gl_viewer_impl
  {
    const GLfloat SplineCurve::SPLINE_SAMPLING_TOLERANCE = 15.0;
    
    Polyline::Polyline(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates)
    {
      _coordinates.resize(2 * coordinates.size());
      
      for(size_t i = 0; i < coordinates.size(); ++i)
      {
        _coordinates[2 * i] = coordinates[i](0);
        _coordinates[2 * i + 1] = coordinates[i](1);
      }
    }
        
    void Polyline::execute(size_t mode)
    {
      glBegin(GL_LINE_STRIP);
      size_t i = 0;
      while(i < _coordinates.size())
      {
        glVertex2f(_coordinates[i], _coordinates[i + 1]);
        i += 2;
      }
      glEnd();
    }
    
    Polygon::Polygon(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates)
    {
      _coordinates.resize(2 * coordinates.size());
      
      for(size_t i = 0; i < coordinates.size(); ++i)
      {
        _coordinates[2 * i] = coordinates[i](0);
        _coordinates[2 * i + 1] = coordinates[i](1);
      }
    }
        
    void Polygon::execute(size_t mode)
    {
      if(_coordinates.size() < 2)
        return;
        
      glBegin(GL_LINE_STRIP);
      size_t i = 0;
      while(i < _coordinates.size())
      {
        glVertex2f(_coordinates[i], _coordinates[i + 1]);
        i += 2;
      }
      glVertex2f(_coordinates[0], _coordinates[1]);
      glEnd();
    }
    
    FillPolygon::FillPolygon(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates)
    {
      _triangulator = 0;
      _coordinates.resize(3 * coordinates.size());
      
      for(size_t i = 0; i < coordinates.size(); ++i)
      {
        _coordinates[3 * i] = coordinates[i](0);
        _coordinates[3 * i + 1] = coordinates[i](1);
        _coordinates[3 * i + 2] = 0.0;
      }
    }
        
    FillPolygon::~FillPolygon()
    {
      if(_triangulator)
        gluDeleteTess(_triangulator);
    }
    
    void FillPolygon::initialize()
    {
      _triangulator = gluNewTess();
      
      gluTessCallback(_triangulator, GLU_TESS_VERTEX, (GLvoid (*) ()) &glVertex3dv);
      gluTessCallback(_triangulator, GLU_TESS_BEGIN, (GLvoid (*) ()) &glBegin);
      gluTessCallback(_triangulator, GLU_TESS_END, (GLvoid (*) ()) &glEnd);
    }
        
    void FillPolygon::execute(size_t mode)
    {
      gluBeginPolygon(_triangulator);
      size_t i = 0;
      while(i < _coordinates.size())
      {
        gluTessVertex(_triangulator, &(_coordinates[i]), &(_coordinates[i]));
        i += 3;
      }
      gluEndPolygon(_triangulator);
    }
    
    void Vertex::execute(size_t mode)
    {
      if(mode == FILE_MODE)
        gl2psPointSize(5.0f);
        
      if(mode == SCREEN_MODE)
        glPointSize(5.0f);
        
      glBegin(GL_POINTS);
      glVertex2f(_coordinates[0], _coordinates[1]);
      glEnd();
    }
    
    Image::Image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> & x_interval, const ublas::fixed_vector<float_t, 2> & y_interval) :
      _image_size(image.size()), _x_interval(x_interval), _y_interval(y_interval), _texture_name(0)
    {
      for(size_t j = 0; j < 2; j++)
      {
        size_t i;
        for (i = 2; i < _image_size(j); i <<= 1);
        _texture_size(j) = i;

        _texture_coordinates(j) = float(_image_size(j)) / float(_texture_size(j));
      }     
      
      _byte_image_data.resize(_texture_size(0) * _texture_size(1) * 3); //array, in das die Bildinformation bertragen wird.

      for(size_t j = 0; j < _image_size(1); j++)
        for(size_t i = 0; i < _image_size(0); i++)
        {
          size_t index = 3*(i + j * _texture_size(0));
          _byte_image_data[index++]= (unsigned char)(image[ublas::fixed_vector<size_t, 2>(i, j)](0));
          _byte_image_data[index++]= (unsigned char)(image[ublas::fixed_vector<size_t, 2>(i, j)](1));
          _byte_image_data[index++]= (unsigned char)(image[ublas::fixed_vector<size_t, 2>(i, j)](2));
        }     
      
      _float_image_data.resize(_image_size(0) * _image_size(1) * 3); //array, in das die Bildinformation bertragen wird.

      for(size_t j = 0; j < _image_size(1); j++)
        for(size_t i = 0; i < _image_size(0); i++)
        {
          size_t index = 3*(i + j * _image_size(0));
          _float_image_data[index++]= img::float_t(image[ublas::fixed_vector<size_t, 2>(i, j)](0));
          _float_image_data[index++]= img::float_t(image[ublas::fixed_vector<size_t, 2>(i, j)](1));
          _float_image_data[index++]= img::float_t(image[ublas::fixed_vector<size_t, 2>(i, j)](2));
        }
    }
    
    void Image::initialize()
    {
      glGenTextures(1, &_texture_name);//gl generiert einen texturnamen
      glBindTexture(GL_TEXTURE_2D, _texture_name);//texturname wird benutzt
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);// lineares downfiltering
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);// linearen upfiltering
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);//die texturkoordinaten sind zwischen 0.0, und 1.0 geclamped
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);//die texturkoordinaten sind zwischen 0.0, und 1.0 geclamped
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, _texture_size(0), _texture_size(1), 0, GL_RGB, GL_UNSIGNED_BYTE, &_byte_image_data[0]);
      
      _byte_image_data.resize(0);
    }
      
    void Image::execute(size_t mode)
    {
      if(mode == SCREEN_MODE)
      {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, _texture_name);
        glBegin(GL_QUADS);

//         glColor3f(1.0f, 1.0f, 1.0f);

        glTexCoord2f(0.0f, 0.0f );
        glVertex2f(_x_interval(0), _y_interval(0));
        
        glTexCoord2f(_texture_coordinates(0), 0.0f);
        glVertex2f(_x_interval(1), _y_interval(0));

        glTexCoord2f(_texture_coordinates(0), _texture_coordinates(1));
        glVertex2f(_x_interval(1), _y_interval(1));

        glTexCoord2f(0.0f, _texture_coordinates(1));
        glVertex2f(_x_interval(0), _y_interval(1));

        glEnd();
        glDisable(GL_TEXTURE_2D);
      }
      
      if(mode == FILE_MODE)
      {
        glRasterPos3f(GLfloat(0.0), GLfloat(0.0), GLfloat(0.0));

        gl2psDrawPixels(_image_size(0), _image_size(1),
                        0, 0, GL_RGB, GL_FLOAT, &_float_image_data[0]);
      }
    }
    
    void DeleteTexture::initialize()
    {
      if(_texture_name)
        glDeleteTextures(1, &_texture_name);
    }
      
    SplineCurve::SplineCurve(const Bspline< ublas::fixed_vector<float_t, 2> > & curve)
    {
      _nurbs_renderer = 0;
      _order = curve.spline_order();
      _coefficients.resize(3 * curve.n_coefficients());
      _knots.resize(curve.n_knots());  
      
      for(size_t i = 0; i < curve.n_coefficients(); ++i)
      {
        _coefficients[3 * i] = curve.coefficient(i)(0);
        _coefficients[3 * i + 1] = curve.coefficient(i)(1);
        _coefficients[3 * i + 2] = 0.0;
      }
          
      for(size_t i = 0; i < curve.n_knots(); ++i)
        _knots[i] = curve.knot(i);
    }
      
    void SplineCurve::initialize()
    {
      _nurbs_renderer = gluNewNurbsRenderer();
      gluNurbsProperty(_nurbs_renderer, GLU_SAMPLING_TOLERANCE, SPLINE_SAMPLING_TOLERANCE);
      gluNurbsProperty(_nurbs_renderer, GLU_DISPLAY_MODE, GLU_FILL);
    }
    
    SplineCurve::~SplineCurve()
    {
      if(_nurbs_renderer)
        gluDeleteNurbsRenderer(_nurbs_renderer);
    }
    
    void SplineCurve::execute(size_t mode)
    {
      gluBeginCurve(_nurbs_renderer);
      gluNurbsCurve(_nurbs_renderer, _knots.size(), &_knots[0], 3, &_coefficients[0], _order, GL_MAP1_VERTEX_3);
      gluEndCurve(_nurbs_renderer);
    }
    
    void LineWidth::execute(size_t mode)
    {
      if(mode == SCREEN_MODE)
        glLineWidth(_line_width);
      
      if(mode == FILE_MODE)
        gl2psLineWidth(_line_width);
    }
    
    SetColor::SetColor(const Color & color)
    {
      for(size_t i = 0; i < 3; ++i)
        _color(i) = (unsigned char)(color(i));
    }
    
    void SetColor::execute(size_t mode)
    {
      glColor3ub(_color(0), _color(1), _color(2));
    }
    
    void OffsetZLayer::execute(size_t mode)
    {
      glTranslatef(GLfloat(0.0), GLfloat(0.0), GLfloat(0.01 * _offset));
    }
    
    void Translate::execute(size_t mode)
    {
      glTranslatef(_translation(0), _translation(1), GLfloat(0.0));
    }
    
    WriteImage::WriteImage(const std::string file_name, size_t width, size_t height, const ublas::fixed_vector<float_t, 2> & x_axis_interval, const ublas::fixed_vector<float_t, 2> & y_axis_interval, int image_format, img::OpenGlViewer & viewer) : 
      _file_name(file_name), _width(width), _height(height), _x_axis_interval(x_axis_interval), _y_axis_interval(y_axis_interval), _viewer(viewer), _executed(false)
    {  
      int local_image_format = image_format;
      
      if(local_image_format == OpenGlViewer::IMAGE_FORMAT_DETERMINE)
        local_image_format = determine_image_format(_file_name);
  
      switch(local_image_format)
      {
      case OpenGlViewer::IMAGE_FORMAT_EPS:
        _image_format = GL2PS_EPS;
        break;
  
      case OpenGlViewer::IMAGE_FORMAT_SVG:
        _image_format = GL2PS_SVG;
        break;
  
      case OpenGlViewer::IMAGE_FORMAT_PDF:
        _image_format = GL2PS_PDF;
        break;
  
      default:
        throw Exception("Exception: Non-supported file format in OpenGlViewer::write_image().");
      }
    }
    
    void WriteImage::execute(size_t mode)
    {
      if(mode == SCREEN_MODE && ! _executed)
      {
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        
        FILE *current_output_file;

        current_output_file = fopen(_file_name.c_str(), "wb");

        if(current_output_file == 0x0)
          throw Exception("Exception: Could not open file " + _file_name + " in WriteImage::execute().");

        int state = GL2PS_OVERFLOW, buffsize = 0;
    
        _viewer.set_viewport(_width, _height, _x_axis_interval, _y_axis_interval);
    
        while(state == GL2PS_OVERFLOW)
        {
          buffsize += 1024*1024;
    
          gl2psBeginPage("Imaging2 OpenGL output", "Imaging2 Class Library", NULL, _image_format,
                        GL2PS_SIMPLE_SORT, GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT,
                        GL_RGBA, 0, NULL, 0, 0, 0, buffsize, current_output_file,
                        _file_name.c_str());
                        
          glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
          glLoadIdentity();        
              
          std::map<size_t, boost::shared_ptr<open_gl_viewer_impl::ObjectInterface> >::iterator object_iter;
          for(object_iter = _viewer._objects.begin(); object_iter != _viewer._objects.end(); ++object_iter)
            object_iter->second->execute(open_gl_viewer_impl::FILE_MODE);
    
          glFlush();
          state = gl2psEndPage();
        }
    
        OpenGlViewer::out.set_viewport(viewport[2], viewport[3]);
    
        fclose(current_output_file);
        
        _executed = true;
      }
    }
    
    int WriteImage::determine_image_format( const std::string & file_name )
    {
      int length = int(file_name.length());
      std::string suffix;
      suffix = file_name.substr(length - 3, length - 1);
  
      if(suffix == "eps")
        return OpenGlViewer::IMAGE_FORMAT_EPS;
  
      if(suffix == "svg")
        return OpenGlViewer::IMAGE_FORMAT_SVG;
  
      if(suffix == "pdf")
        return OpenGlViewer::IMAGE_FORMAT_PDF;
  
      return OpenGlViewer::IMAGE_FORMAT_CAN_NOT_DETERMINE;
    }
  }
}

