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

#ifndef LAPACK_LINEARALGEBRA_H
#define LAPACK_LINEARALGEBRA_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \ingroup la
      <tt>\#include <lapack/linear_algebra.hpp></tt>
  
      Computes the eigenvectors and the corresponding eigenvalues of a symmetric and positive semidefinite matrix \em A. The \em eigenvalues are sorted from the largest to the smallest. The columns of \em eigenvectors are ordered in the same way.
      If \f$\textrm{diag}(\sigma_1, \ldots, \sigma_n)\f$ is the matrix with components of \em eigenvalues in the diagonal and \f$V\f$ denotes the matrix \em eigenvectors, then 
      \f[
        A V = V \textrm{diag}(\sigma_1, \ldots, \sigma_n)
      \f]
      holds.
  */
  void eigensystem(const ublas::matrix<float_t> & A, ublas::matrix<float_t> & eigenvectors, ublas::vector<float_t> & eigenvalues);
  
  /** \ingroup la
      <tt>\#include <lapack/linear_algebra.hpp></tt>
      
      For a symmetric and positive semidefinite matrix \em A, this function computes \em root such that if denote \em root as \f$R\f$ the relation
      \f[
        A = RR^T 
      \f]
      holds.
  */
  void symmetric_square_root(const ublas::matrix<float_t> & A, ublas::matrix<float_t> & root);
  
  /** \ingroup la
      <tt>\#include <lapack/linear_algebra.hpp></tt>
      
      For a symmetric and positive semidefinite matrix \em A, this function computes \em root such that if denote \em root as \f$R\f$ the relation
      \f[
        A = RR 
      \f]
      holds.
  */
  void square_root(const ublas::matrix<float_t> & A, ublas::matrix<float_t> & root);  
  
  /** \ingroup la
      <tt>\#include <lapack/linear_algebra.hpp></tt>
      
      For a symmetric and positive semidefinite matrix \em A, this function computes \em root such that if denote \em root as \f$R\f$ the relation
      \f[
        A = R^{-1}R^{-1} 
      \f]
      holds.
  */
  void inverse_square_root(const ublas::matrix<float_t> & A, ublas::matrix<float_t> & root);
  
  /** \ingroup la
      <tt>\#include <lapack/linear_algebra.hpp></tt>
      
      Computes the \em inverse of a non-singular matrix \em A.
  */
  void inverse(const ublas::matrix<float_t> & A, ublas::matrix<float_t> & inverse);
}

#endif
