// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_CIO_H
#define CORE_CIO_H

#include <core/imaging2.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <set>
#include <map>


namespace imaging
{
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  template <class element_t>
  std::ostream & operator<<(std::ostream & out, const std::vector<element_t> & vector)
  {
    out << "[" << vector.size() << "](";
    for(size_t i = 0; i < vector.size() - 1; ++i)
      out << vector[i] << ",";
      
    out << vector[vector.size() - 1] << ")";
      
    return out;
  }
  
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  template <class key_t, class data_t>
  std::ostream & operator<<(std::ostream & out, const std::map<key_t, data_t> & map)
  {
    typename std::map<key_t, data_t>::const_iterator iter;
    
    for(iter = map.begin(); iter != map.end(); ++iter)
      out << "(" << iter->first << "," << iter->second << ") ";
      
    return out;
  } 
  
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  template <class data_t>
  std::ostream & operator<<(std::ostream & out, const std::set<data_t> & set)
  {
    typename std::set<data_t>::const_iterator iter;
    
    for(iter = set.begin(); iter != set.end(); ++iter)
      out << *iter << ", ";
      
    return out;
  } 
  
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  template <class key_t, class data_t>
  std::ostream & operator<<(std::ostream & out, const std::pair<key_t, data_t> & pair)
  {
    out << "(" << pair.first << ", " << pair.second << ")";
      
    return out;
  } 
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  void output_matrix_matlab_style(std::ostream & out, const ublas::matrix<float_t> & matrix);
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  void output_vector_matlab_style(std::ostream & out, const ublas::vector<float_t> & vector);
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  void output_matrix_gnuplot_style(std::ostream & out, const ublas::matrix<float_t> & matrix);
  
  /** \ingroup core
      <tt>\#include <core/cio.hpp></tt>
  */
  void output_vector_gnuplut_style(std::ostream & out, const ublas::vector<float_t> & vector);
}


#endif
