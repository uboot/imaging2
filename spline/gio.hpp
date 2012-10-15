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

#ifndef SPLINE_GIO_H
#define SPLINE_GIO_H

#include <spline/Bspline.hpp>
#include <spline/PeriodicBspline.hpp>
#include <graphics/GraphicsInterface.hpp>


namespace imaging
{
  /** \ingroup spline 
      <tt>\#include <spline/gio.hpp></tt>
    
      Draws the 2-dimensional B-Spline curve \em curve in the graphics output. */
  GraphicsInterface & operator<<(GraphicsInterface & out, const Bspline< ublas::fixed_vector<float_t, 2> > & curve);
  
  /** \ingroup spline 
      <tt>\#include <spline/gio.hpp></tt>
    
      Draws the 2-dimensional B-Spline curve \em curve in the graphics output. */
  GraphicsInterface & operator<<(GraphicsInterface & out, const PeriodicBspline< ublas::fixed_vector<float_t, 2> > & curve);
}


#endif
