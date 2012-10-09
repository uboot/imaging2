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

#ifndef CORE_MATRIXUTILITIES_H
#define CORE_MATRIXUTILITIES_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      Returns the determinant of \em A.
  */
  float_t determinant(const ublas::fixed_matrix<float_t, 2, 2> & A);
  
  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      Returns the inverse of \em A.
  */
  ublas::fixed_matrix<float_t, 2, 2> inverse(const ublas::fixed_matrix<float_t, 2, 2> & A);
  
  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      Returns the determinant of \em A.
  */
  float_t determinant(const ublas::fixed_matrix<float_t, 3, 3> & A);
  
  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      Returns the inverse of \em A.
  */
  ublas::fixed_matrix<float_t, 3, 3> inverse(const ublas::fixed_matrix<float_t, 3, 3> & A);

  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      Returns the inverse of \em A.
  */
  ublas::fixed_matrix<float_t, 1, 1> inverse(const ublas::fixed_matrix<float_t, 1, 1> & A);
  
  
  /** \ingroup core
      <tt>\#include <core/matrix_utilities.hpp></tt>
      
      For the angle \em alpha, this function returns the rotation matrix
      \f[
      \left(
      \begin{array}{cc}
      \cos \alpha & -\sin \alpha \\
      \sin \alpha & \cos \alpha
      \end{array}
      \right)\,.
      \f]
  */
  ublas::fixed_matrix<float_t, 2, 2> rotation_matrix(float_t alpha);
}

#endif
