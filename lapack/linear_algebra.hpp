// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
