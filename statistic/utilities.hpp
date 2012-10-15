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
