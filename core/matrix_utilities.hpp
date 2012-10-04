// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
