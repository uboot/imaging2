// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
