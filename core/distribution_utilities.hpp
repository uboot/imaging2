// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_DISTRIBUTIONUTILITIES_H
#define CORE_DISTRIBUTIONUTILITIES_H

#include <core/imaging2.hpp>
namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/distribution_utilities.hpp></tt>
      
      Returns a sample from the uniform distribution on the interval [0, 1].
  */
  float_t uniform_distribution();
  
  /** \ingroup core
      <tt>\#include <core/distribution_utilities.hpp></tt>
      
      Returns a sample from the uniform distribution on the interval [-1, 1].
  */
  float_t symmetric_uniform_distribution();
  
  /** \ingroup core
      <tt>\#include <core/distribution_utilities.hpp></tt>
      
      Returns a sample from the <em>N(0, 1)</em>-normal distribution.
  */
  float_t normal_distribution();
}

#endif
