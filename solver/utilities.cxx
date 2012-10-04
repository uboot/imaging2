#include <solver/utilities.hpp>

namespace imaging
{
  /** \cond */
  void sparse2raw(const ublas::compressed_matrix<float_t> & matrix, std::vector<float_t> & values, std::vector<int> & column_indices, std::vector<int> & block_indices)
  {
    ublas::compressed_matrix<float_t>::const_iterator1 iter1;
    ublas::compressed_matrix<float_t>::const_iterator2 iter2;

    std::size_t n_elements = 0;

    for(iter1 = matrix.begin1(); iter1 != matrix.end1(); ++iter1)
      for(iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        n_elements++;

    values.resize(n_elements);
    column_indices.resize(n_elements);
    block_indices.resize(matrix.size1() + 1);

    block_indices[0] = 1;

    std::size_t index = 0;

    for(iter1 = matrix.begin1(); iter1 != matrix.end1(); ++iter1)
    {
      for(iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
      {
        values[index] = *iter2;
        column_indices[index] = iter2.index2() + 1;
        index++;
      }
      block_indices[iter1.index1() + 1] = index + 1;
    }
  }
  /** \endcond */
}

