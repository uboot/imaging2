// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CGSOLVER_H
#define CGSOLVER_H

#include <solver/SolverInterface.hpp>

namespace imaging
{

  /** \ingroup solver
      \brief Wrapper class for the CG solver of provided by ITPACK.
      
      This class implements a conjugated gradients solver by interfacing <a href="http://rene.ma.utexas.edu/CNA/ITPACK">ITPACK</a>.
  */
  class CgSolver : public SolverInterface
  {
    size_t _n_max_iterations;

    void report_error(int error_flag) const;

  public:
    /** Constructs a CG solver. The solver iterates at most \em n_max_iterations times. */
    CgSolver(size_t n_max_iterations = 10000) : _n_max_iterations(n_max_iterations) {}

    void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const;
  };  
}

#endif
