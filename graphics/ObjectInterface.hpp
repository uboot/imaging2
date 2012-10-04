// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef GRAPHICS_OPENGLVIEWERIMPL_OBJECTINTERFACE_H
#define GRAPHICS_OPENGLVIEWERIMPL_OBJECTINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  namespace open_gl_viewer_impl
  {
    enum draw_modes { SCREEN_MODE, FILE_MODE };
    
    class ObjectInterface
    {
    public:
      virtual ~ObjectInterface() {}
      
    /* GLUT thread only BEGIN */
      virtual void initialize() {}
      virtual void execute(size_t mode) = 0;
    /* GLUT thread only END */
    };
  }
}

#endif
