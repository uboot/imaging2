// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
