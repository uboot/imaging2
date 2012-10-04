// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_SHAPEENERGYINTERFACE_H
#define SHAPE_SHAPEENERGYINTERFACE_H

#include <minimize/EnergyInterface.hpp>


namespace imaging
{
  /** \ingroup shape 
      \brief Abstract class interface for energies depending on vector valued input data which define an input shape.
      
      In addition to the requirements of EnergyInterface, classes implementing ShapeEnergyInterface provide a current shape computed from the current argument. Minimizers can access this shape to obtain geometric information (via BoundaryDiscretizer) about the current argument. 
  */
  template <class shape_t>
  class ShapeEnergyInterface : public EnergyInterface
  {
  public:
  
    /** Returns the current shape. This function should \em not actually compute the current shape but return the cached result of the last call to set_argument()! */ 
    virtual const shape_t & current_shape() const = 0;
  };
}

#endif
