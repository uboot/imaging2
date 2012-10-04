// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_FUNCTIONINTERFACE_H
#define CORE_FUNCTIONINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class template for functions.
       
      Objects implementing this class template map objects of type \em A to objects of type \em B.
  */
  template <class A, class V>
  class FunctionInterface
  {
  public:
    virtual ~FunctionInterface() {}
  
    /** Evaluates the function at \em argument. */
    virtual void evaluate(const A & argument, V & value) = 0;
  };
}

#endif
