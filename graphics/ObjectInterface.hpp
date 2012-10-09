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
