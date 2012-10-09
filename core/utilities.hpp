/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CORE_UTILITIES
#define CORE_UTILITIES


#include <iostream>
#include <sstream>
#include <math.h>
#include <core/imaging2.hpp>

#include <boost/lexical_cast.hpp>


namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Defines Pi. */
  const double PI = 3.1415926;
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Defines the square root of 2. */
  const double SQUARE_ROOT_2 = 1.414213562373095049;
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Defines the square root of 2. */
  const double SQUARE_ROOT_3 = 1.73205080757;

  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns max(a, b). */
  template <class num_t>
  inline num_t max(const num_t a, const num_t b) { return a > b ? a : b; }

  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns min(a, b). */
  template <class num_t>
  inline num_t min(const num_t a, const num_t b) { return a < b ? a : b; }

  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns a<sup>2</sup>. */
  template <class num_t>
  inline num_t square(num_t a) { return a * a; }

  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns a<sup>b</sup>, where b is a non-negative integer. */
  template <class num_t>
  inline num_t power(num_t a, size_t b)
  { 
    if(b == 0)
      return num_t(1.0);
    else
      return power(a, b - 1) * a;
  }
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns the absolute value of \em a. */
  template <class num_t>
  inline num_t abs(num_t a)
  {
    return a < 0 ? -a : a;
  }
  
  template <>
  inline float abs(float a)
  {
    return fabs(a);
  }
  
  template <>
  inline double abs(double a)
  {
    return fabs(a);
  }
  
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns the sign of \em a. */
  template <class num_t>
  inline num_t sign(num_t a)
  {
    if(a == num_t(0))
      return num_t(0);
      
    return a < 0 ? num_t(-1) : num_t(1);
  }
  
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
   
      Returns 1.0 if \em base equals \em t and 0.0 otherwise. */
  float_t delta(size_t base, size_t t);
  
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns the clockwise difference from \em angle_1 to \em angle_2 in the interval [0, 2 Pi]. */
  float_t clockwise_difference(float_t angle_1, float_t angle_2);
  
  
  /** \ingroup core
      <tt>\#include <core/utilities.hpp></tt>
      
      Returns the counterclockwise difference from \em angle_1 to \em angle_2 in the interval [0, 2 Pi]. */
  float_t counter_clockwise_difference(float_t angle_1, float_t angle_2);
}

#endif
