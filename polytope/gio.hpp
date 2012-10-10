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

#ifndef POLYTOPE_GIO_H
#define POLYTOPE_GIO_H

#include <polytope/Polygon.hpp>
#include <graphics/GraphicsInterface.hpp>


namespace imaging
{ 
  /** \ingroup polytope 
      <tt>\#include <polytope/gio.hpp></tt> */
  GraphicsInterface & operator<<(GraphicsInterface & out, const Polygon & polygon);
  
  /** \ingroup polytope 
      <tt>\#include <polytope/gio.hpp></tt> */
  GraphicsInterface & operator<<(GraphicsInterface & out, const SimplePolygon & polygon);
}


#endif
