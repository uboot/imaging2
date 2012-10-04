// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SOLVER_UTILTIES_H
#define SOLVER_UTILTIES_H

#include <core/imaging2.hpp>

namespace imaging
{
  /** \cond */
  void sparse2raw(const ublas::compressed_matrix<float_t> & matrix, std::vector<float_t> & values, std::vector<int> & column_indices, std::vector<int> & block_indices);
  /** \endcond */
}

#endif
