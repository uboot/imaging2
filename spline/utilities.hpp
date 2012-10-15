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

#ifndef SPLINE_UTILITIES_H
#define SPLINE_UTILITIES_H
 
#include <spline/Bspline.hpp>
#include <spline/PeriodicBspline.hpp>

namespace imaging
{ 
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Sets the <em>i</em>-th spline coefficient of \em spline to 1 and all the others to 0.
  */
  void basis_spline(size_t i, Bspline<float_t> & spline);
  
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Sets the <em>i</em>-th spline coefficient of \em spline to 1 and all the others to 0. */
  void basis_spline(size_t i, PeriodicBspline<float_t> & spline);
  
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Computes the matrix which is the result of the spline interpolation problem. Assume a B-spline function
      \f[
        s(\tau) = \sum_{i = 0}^n a_i b_i(\tau)
      \f]
      and consider the interpolation problem
      \f[
        s(\tau_k) = y_k\,,\quad 0 \leq k \leq N\,.
      \f]
      Then the least squares solution of this problem is defined by the system of linear equations
      \f[
        \sum_{i=0}^n a_i \underbrace{\sum_{k=0}^N b_i(\tau_k)b_j(\tau_k)}_{A_{ij}} =
        \sum_{j=0}^N \underbrace{b_j(\tau_k)}_{B_{jk}} y_k\,,\quad 0 \leq j \leq n\,.
      \f]
      If \f$s\f$ is given by \em spline and \f$(\tau_k)_{0\leq k \leq N}\f$ by \em values, then this function writes \f$A\f$ to \em matrix.
  */
  void interpolation_matrix(const Bspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
  
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Computes the matrix which is the result of the spline interpolation problem. Assume a B-spline function
      \f[
        s(\tau) = \sum_{i = 0}^n a_i b_i(\tau)
      \f]
      and consider the interpolation problem
      \f[
        s(\tau_k) = y_k\,,\quad 0 \leq k \leq N\,.
      \f]
      Then the least squares solution of this problem is defined by the system of linear equations
      \f[
        \sum_{i=0}^n a_i \underbrace{\sum_{k=0}^N b_i(\tau_k)b_j(\tau_k)}_{A_{ij}} =
        \sum_{j=0}^N \underbrace{b_j(\tau_k)}_{B_{jk}} y_k\,,\quad 0 \leq j \leq n\,.
      \f]
      If \f$s\f$ is given by \em spline and \f$(\tau_k)_{0\leq k \leq N}\f$ by \em values, then this function writes \f$A\f$ to \em matrix.
  */
  void interpolation_matrix(const PeriodicBspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
  
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Computes the matrix which is the result of the spline interpolation problem. Assume a B-spline function
      \f[
        s(\tau) = \sum_{i = 0}^n a_i b_i(\tau)
      \f]
      and consider the interpolation problem
      \f[
        s(\tau_k) = y_k\,,\quad 0 \leq k \leq N\,.
      \f]
      Then the least squares solution of this problem is defined by the system of linear equations
      \f[
        \sum_{i=0}^n a_i \underbrace{\sum_{k=0}^N b_i(\tau_k)b_j(\tau_k)}_{A_{ij}} =
        \sum_{j=0}^N \underbrace{b_j(\tau_k)}_{B_{jk}} y_k\,,\quad 0 \leq j \leq n\,.
      \f]
      If \f$s\f$ is given by \em spline and \f$(\tau_k)_{0\leq k \leq N}\f$ by \em values, then this function writes \f$B\f$ to \em matrix.
  */
  void basis_spline_matrix(const PeriodicBspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
  
  /** \ingroup spline
      <tt>\#include <spline/utilities.hpp></tt>
     
      Computes the matrix which is the result of the spline interpolation problem. Assume a B-spline function
      \f[
        s(\tau) = \sum_{i = 0}^n a_i b_i(\tau)
      \f]
      and consider the interpolation problem
      \f[
        s(\tau_k) = y_k\,,\quad 0 \leq k \leq N\,.
      \f]
      Then the least squares solution of this problem is defined by the system of linear equations
      \f[
        \sum_{i=0}^n a_i \underbrace{\sum_{k=0}^N b_i(\tau_k)b_j(\tau_k)}_{A_{ij}} =
        \sum_{j=0}^N \underbrace{b_j(\tau_k)}_{B_{jk}} y_k\,,\quad 0 \leq j \leq n\,.
      \f]
      If \f$s\f$ is given by \em spline and \f$(\tau_k)_{0\leq k \leq N}\f$ by \em values, then this function writes \f$B\f$ to \em matrix.
  */
  void basis_spline_matrix(const Bspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
  
  /** \cond */
  namespace spline_utilities_impl
  { 
    template <class spline_t>
    void basis_spline(size_t i, spline_t & spline);
    
    template <class spline_t>
    void interpolation_matrix(const spline_t & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
    
    template <class spline_t>
    void basis_spline_matrix(const spline_t & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix);
  }
  /** \endcond */
}
 
 
#endif
