// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_VECTORUTILITIES_H
#define CORE_VECTORUTILITIES_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/vector_utilities.hpp></tt>
      
      Converts the polar coordinates <em>(radius, angle)</em> to cartesian coordinates and returns the result. Mathematically:
      \f[
      (r, \alpha) \mapsto (r \cos\alpha, r\sin\alpha)
      \f]
  */
  const ublas::fixed_vector<float_t, 2> polar2cartesian(float_t radius, float_t angle);
  
  
  /** \ingroup core
      <tt>\#include <core/vector_utilities.hpp></tt>
      
      Returns the length of \em v.
  */
  float_t radius(const ublas::fixed_vector<float_t, 2> & v);
  
  
  /** \ingroup core
      <tt>\#include <core/vector_utilities.hpp></tt>
      
      Returns the polar angle of \em v.
  */
  float_t angle(const ublas::fixed_vector<float_t, 2> & v);
}

#endif
