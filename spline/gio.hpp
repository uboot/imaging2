// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
