// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SOLVER_SOLVERINTERFACE_H
#define SOLVER_SOLVERINTERFACE_H

#include <core/imaging2.hpp>
                    
namespace imaging
{                                      
  /** \ingroup solver
      \brief Abstract class interface for solver for sparse systems of linear equations.
  */
  class SolverInterface
  {

  public:
    virtual ~SolverInterface() {}
  
    /** Solves the system of equations defined by the matrix \em eqs and the vector \em rhs and writes the solution to \em result. The vector \em result is automatically resized to dimension of the system. If the sizes of the input data are such that the system cannot be solved an Exception is thrown. */
    virtual void solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const = 0;
  };

}

#endif
