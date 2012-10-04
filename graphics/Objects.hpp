// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef GRAPHICS_OPENGLVIEWERIMPL_OBJECTS_H
#define GRAPHICS_OPENGLVIEWERIMPL_OBJECTS_H

#include <graphics/ObjectInterface.hpp>
#include <image/Color.hpp>
#include <spline/Bspline.hpp>
#include <glut.h>

namespace imaging
{
  class OpenGlViewer;
  
  namespace open_gl_viewer_impl
  {
    class Polyline : public ObjectInterface
    {
      std::vector<GLfloat> _coordinates;
      
    public:
      Polyline(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates);
      void execute(size_t mode);
    };
    
    class Polygon : public ObjectInterface
    {
      std::vector<GLfloat> _coordinates;
      
    public:
      Polygon(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates);
      void execute(size_t mode);
    };
    
    class FillPolygon : public ObjectInterface
    {
      std::vector<GLdouble> _coordinates;
      GLUtriangulatorObj * _triangulator;
      
    public:
      FillPolygon(const std::vector< ublas::fixed_vector<float_t, 2> > coordinates);
      ~FillPolygon();
      void initialize();
      void execute(size_t mode);
    };
    
    class Vertex : public ObjectInterface
    {
      ublas::fixed_vector<GLfloat, 2> _coordinates;
      
    public:
      Vertex(const ublas::fixed_vector<float_t, 2> & vertex) :
        _coordinates(vertex) {}
      void execute(size_t mode);
    };

    class Image : public ObjectInterface
    {
      ublas::fixed_vector<GLint, 2> _image_size;
      ublas::fixed_vector<GLfloat, 2> _texture_coordinates;
      ublas::fixed_vector<GLfloat, 2> _x_interval;
      ublas::fixed_vector<GLfloat, 2> _y_interval;
      
      ublas::fixed_vector<GLuint, 2> _texture_size;
      GLuint _texture_name;
      std::vector<unsigned char> _byte_image_data;
      std::vector<GLfloat> _float_image_data;
      
    public:
      Image(const ColorImage2d & image, const ublas::fixed_vector<float_t, 2> & x_interval, const ublas::fixed_vector<float_t, 2> & y_interval);
      GLuint texture_name() const { return _texture_name; }
      void initialize();
      void execute(size_t mode);
    };
    
    class DeleteTexture : public ObjectInterface
    {
      GLuint _texture_name;
    
    public:
      DeleteTexture(GLuint texture_name) : _texture_name(texture_name) {}
      void initialize();
      void execute(size_t mode) {}
    };
    
    class ClearPipeline : public ObjectInterface
    { 
    public:
      void execute(size_t mode) {}
    };

    class SplineCurve : public ObjectInterface
    {
      static const GLfloat SPLINE_SAMPLING_TOLERANCE;
      std::vector<GLfloat> _knots;
      std::vector<GLfloat> _coefficients;
      GLint _order;
      GLUnurbs* _nurbs_renderer;
    
    public:  
      SplineCurve(const Bspline< ublas::fixed_vector<float_t, 2> > & spline_curve);
      ~SplineCurve();
      void initialize();
      void execute(size_t mode);
    };
    
    class LineWidth : public ObjectInterface
    {
      GLfloat _line_width;
      
    public:
      LineWidth(float_t line_width) : _line_width(line_width) {}
      void execute(size_t mode);
    };
    
    class SetColor : public ObjectInterface
    {
      ublas::fixed_vector<GLubyte, 3> _color;
      
    public:
      SetColor(const Color & color);
      void execute(size_t mode);
    };
    
    class OffsetZLayer : public ObjectInterface
    {
      GLint _offset;
      
    public:
      OffsetZLayer(int offset) : _offset(offset) {}
      void execute(size_t mode);
    };
    
    class Translate : public ObjectInterface
    {
      ublas::fixed_vector<GLfloat, 2> _translation;
      
    public:
      Translate(const ublas::fixed_vector<float_t, 2> & translation) :_translation(translation) {}
      void execute(size_t mode);
    };
    
    class WriteImage : public ObjectInterface
    {
      std::string _file_name;
      int _width, _height;
      GLint _image_format;
      img::OpenGlViewer & _viewer;
      bool _executed;
      const ublas::fixed_vector<float_t, 2> _x_axis_interval;
      const ublas::fixed_vector<float_t, 2> _y_axis_interval;
      
      
      static int determine_image_format( const std::string & file_name );
      
    public:
      WriteImage(const std::string file_name, size_t width, size_t height, const ublas::fixed_vector<float_t, 2> & x_axis_interval, const ublas::fixed_vector<float_t, 2> & y_axis_interval, int image_format, img::OpenGlViewer & viewer);
      
      void execute(size_t mode);
    };
      
  }
}

#endif
