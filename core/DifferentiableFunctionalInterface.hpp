// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_DIFFERENTIABLEFUNCTIONALINTERFACE_H
#define CORE_DIFFERENTIABLEFUNCTIONALINTERFACE_H

#include <core/FunctionalInterface.hpp>

namespace imaging
{
  /** \ingroup core 
      \brief Abstract class interface of differentiable functionals depending on vector valued input data. */
  class DifferentiableFunctionalInterface : public FunctionalInterface
  {
  public:
    using FunctionalInterface::operator();
    
    virtual float_t operator()(const ublas::vector<float_t> & x, ublas::vector<float_t> & gradient) = 0;
  };
}

#endif
