// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef STATISTIC_UTILITIES_H
#define STATISTIC_UTILITIES_H

#include <core/imaging2.hpp>


namespace imaging
{
  /** \ingroup statistic
      <tt>\#include <statistic/utilities.hpp></tt>
      
      Computes the mean of each column of \em data and writes it to \em result. 
      The vector \em result is resized to the number of columns of \em data upon return.
  */
  void mean(const ublas::matrix<float_t> & data, ublas::vector<float_t> & result);
  
  /** \ingroup statistic
      <tt>\#include <statistic/utilities.hpp></tt>
      
      Computes the variances of the columns of \em data and writes them to \em result. 
      The vector \em result is resized to the number of columns of \em data upon return.
  */
  void var(const ublas::matrix<float_t> & data, ublas::vector<float_t> & result);
  
  /** \ingroup statistic
      <tt>\#include <statistic/utilities.hpp></tt>
      
      Returns the mean of \em data.
  */
  float_t mean(const ublas::vector<float_t> & data);
  
  /** \ingroup statistic
      <tt>\#include <statistic/utilities.hpp></tt>
      
      Returns the variance \em data.
  */
  float_t var(const ublas::vector<float_t> & data);
}

#endif
