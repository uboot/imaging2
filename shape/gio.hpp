// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_GIO_H
#define SHAPE_GIO_H

#include <shape/Circle.hpp>
#include <shape/BsplineShape.hpp>
#include <graphics/GraphicsInterface.hpp>


namespace imaging
{   
  /** \ingroup shape
      <tt>\#include <shape/gio.hpp></tt> */
  GraphicsInterface & operator<<(GraphicsInterface & out, const Circle & circle);
  
  /** \ingroup shape
      <tt>\#include <shape/gio.hpp></tt> */
  GraphicsInterface & operator<<(GraphicsInterface & out, const BsplineShape & spline_shape);
}


#endif
