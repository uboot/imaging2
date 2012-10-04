// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SOLVER_LUSOLVER_H
#define SOLVER_LUSOLVER_H

#include <solver/SolverInterface.hpp>
                    
namespace imaging
{                                      
  /** \ingroup solver
      \brief Wrapper class for the LU solver provided by SPOOLES.
      
      This class implements a conjugated gradients solver by interfacing <a href="http://www.netlib.org/linalg/spooles/spooles.2.2.html">SPOOLES</a>.
  */
  class LuSolver : public SolverInterface
  {

  public:
    LuSolver() {}
    
    void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const;
  };

}

#endif
