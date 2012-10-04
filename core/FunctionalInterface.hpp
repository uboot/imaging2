// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_FUNCTIONALINTERFACE_H
#define CORE_FUNCTIONALINTERFACE_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class interface for functionals depending on vector valued input data.
  */
  class FunctionalInterface
  {
  public:
    virtual ~FunctionalInterface() {}
  
    /** Evaluates the functional at \em x. */
    virtual float_t operator()(const ublas::vector<float_t> & x) = 0;
    
    /** Returns the dimension of the input argument for this functional. In other words, this member returns the dimension of the space the functional is defined on. */
    virtual std::size_t dimension() const = 0;
  };
}

#endif
