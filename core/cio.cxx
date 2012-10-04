#include <core/cio.hpp>


namespace imaging
{
  void output_matrix_matlab_style(std::ostream & out, const ublas::matrix<float_t> & matrix)
  {
    out << "[";
    for(ublas::matrix<float_t>::const_iterator1 iter1 = matrix.begin1(); iter1 != matrix.end1(); ++iter1)
    {
      for(ublas::matrix<float_t>::const_iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        out << *iter2 << " ";
      
      out << "; ";
    }
    out << "]";
  } 
  
  
  void output_vector_matlab_style(std::ostream & out, const ublas::vector<float_t> & vector)
  {
    out << "[";
    for(ublas::vector<float_t>::const_iterator iter = vector.begin(); iter != vector.end(); ++iter)
      out << *iter << " ";
    out << "]";
  }
 
  void output_matrix_gnuplot_style(std::ostream & out, const ublas::matrix<float_t> & matrix)
  {
    for(ublas::matrix<float_t>::const_iterator1 iter1 = matrix.begin1(); iter1 != matrix.end1(); ++iter1)
    {
      for(ublas::matrix<float_t>::const_iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        out << *iter2 << " ";
      
      out << "\n";
    }
  } 
  
  void output_vector_gnuplut_style(std::ostream & out, const ublas::vector<float_t> & vector)
  {
    for(ublas::vector<float_t>::const_iterator iter = vector.begin(); iter != vector.end(); ++iter)
      out << *iter << " ";
  }
  
}

