#include <fem/Assembler.hpp>

namespace imaging
{
  void Assembler::clear_matrix(ublas::compressed_matrix<float_t> & matrix)
  {
    ublas::compressed_matrix<float_t>::iterator1 iter1;
    ublas::compressed_matrix<float_t>::iterator2 iter2;
    
    for(iter1 = matrix.begin1(); iter1 != matrix.end1(); ++iter1)
      for(iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = 0;
  }
}

