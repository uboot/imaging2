#include <solver/LuSolver.hpp>

#include <solver/utilities.hpp>

extern "C"
{
  int spooles_lu_solve(double const *values, int const *column_indices,
                       int n_block_indices, int const *block_indices,
                       double const *rhs, double *result);
}
                     
namespace imaging
{
  void LuSolver::solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const
  {
  
  if(eqs.size1() != eqs.size2() || eqs.size2() != rhs.size())
    throw Exception("Exception: Dimensions do not agree in LuSolver::solve().");

  std::vector<float_t> values;
  std::vector<int> column_indices;
  std::vector<int> block_indices;

  sparse2raw(eqs, values, column_indices, block_indices);

  result.resize(rhs.size());
  
  spooles_lu_solve(&values[0], &column_indices[0], 
                   block_indices.size(), &block_indices[0],
                   &rhs[0], &result[0]);
  }
}

