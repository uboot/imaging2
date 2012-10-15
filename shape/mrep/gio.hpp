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

#ifndef SHAPE_MREP_GIO_H
#define SHAPE_MREP_GIO_H

#include <graphics/GraphicsInterface.hpp>
#include <shape/mrep/MrepModel2d.hpp>
#include <shape/mrep/PolygonModel2d.hpp>


namespace imaging
{  
  /** \ingroup mrep
      <tt>\#include <shape/mrep/gio.hpp></tt>
  */
  GraphicsInterface & operator<<(GraphicsInterface & out, const MrepModel2d & model);
  
  
  /** \ingroup mrep
      <tt>\#include <shape/mrep/gio.hpp></tt>
  */
  GraphicsInterface & operator<<(GraphicsInterface & out, const PolygonModel2d & model);
}


#endif
