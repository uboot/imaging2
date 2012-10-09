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
