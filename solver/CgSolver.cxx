#include <solver/CgSolver.hpp>
#include <solver/utilities.hpp>
#include <core/MessageInterface.hpp>

#include <typeinfo>

extern "C" {
void jcg_ (int *NN, int *IA, int *JA, double *A, double *RHS, double *U, 
	   int *IWKSP, int *NW, double *WKSP, int *IPARM, double *RPARM, int *IERR);

void dfault_ (int *IPARM, double *RPARM);
}


namespace imaging
{
  void CgSolver::solve(const ublas::compressed_matrix<float_t> & eqs, const ublas::vector<float_t> & rhs, ublas::vector<float_t> & result) const
  {
    if(eqs.size1() != eqs.size2() || eqs.size2() != rhs.size())
      throw Exception("Exception: Dimensions do not agree in CgSolver::solve().");

    std::vector<float_t> values;
    std::vector<int> column_indices;
    std::vector<int> block_indices;

    sparse2raw(eqs, values, column_indices, block_indices);

    result.resize(rhs.size());

    int NN, *IA, *JA, NW, NW_bkp, *IWKSP = 0, IERR, IPARM[12];
    float_t *ptr_A, *U, *WKSP = 0, RPARM[12];
    float_t const * RHS;

    NN = eqs.size1();

    IA = & block_indices[0];
    JA = & column_indices[0];
    ptr_A = & values[0];
    RHS = & rhs[0];
    U = & result[0];


    if (typeid(float_t) == typeid(double))
      dfault_ (IPARM, (double*)RPARM);
    else throw Exception
      ("Exception: Cannot solve sparse linear system for other type than 'double' in CgSolver::solve()!");

    IPARM[0] = _n_max_iterations;
    IPARM[4] = true;
    NW =  NN * 4 + 4 * IPARM[0];
    NW_bkp = NW;

    IWKSP = new int[NN * 3];
    WKSP = new float_t[NW];

    if (IWKSP && WKSP)
    {
      if (typeid(float_t) == typeid(double))
        jcg_ (&NN, IA, JA, (double*)ptr_A, (double*)RHS, (double*)U,
              IWKSP, &NW, (double*)WKSP, IPARM, (double*)RPARM, &IERR);
      else throw Exception
      ("Exception: Cannot solve sparse linear system for other type than 'double' in CgSolver::solve()!");

      if (IWKSP) delete[] IWKSP;

      if (NW <= NW_bkp)
        delete[] WKSP;
      else
        throw Exception
        ("Exception: Size of WKSP increased during execution in CgSolver::solve()!");

      report_error(IERR);
    }
    else
      throw Exception
        ("Exception: Could not allocate workspace in CgSolver::solve()!");
  }


  void CgSolver::report_error(int error_flag) const
  {
    switch(error_flag)
    {
//     case 0:
//       MessageInterface::out("CgSolver: No errors occured.", MessageInterface::DEBUG_ONLY);
//       break;
    case 11:
      throw Exception("CgSolver: Invalid order of the system!!");
    case 12:
      throw Exception("CgSolver: Workspace array WKSP(*) is not large enough!!");
    case 13:
      MessageInterface::out("CgSolver (Warning): Failure to converge in  IPARM(0) iterations!!", MessageInterface::DEBUG_ONLY);
      break;
    case 14:
      throw Exception("CgSolver: Invalid order of the black subsystem !!");
    case 101:
    case 401:
      MessageInterface::out("CgSolver (Warning): A diagonal element is not positive!!", MessageInterface::DEBUG_ONLY);
      break;
    case 102:
    case 402:
      MessageInterface::out("CgSolver (Warning): No diagonal element in a row!!", MessageInterface::DEBUG_ONLY);
      break;
    case 201:
      MessageInterface::out("CgSolver (Warning): Red-black indexing is not possible!!", MessageInterface::DEBUG_ONLY);
      break;
    case 301:
      MessageInterface::out("CgSolver (Warning): No entry in a row of the original matrix!!", MessageInterface::DEBUG_ONLY);
      break;
    case 302:
      MessageInterface::out("CgSolver (Warning): No entry in a row of the permuted matrix!!", MessageInterface::DEBUG_ONLY);
      break;
    case 303:
      MessageInterface::out("CgSolver (Warning): Sorting error in a row of the permuted matrix!!", MessageInterface::DEBUG_ONLY);
      break;
    case 501:
      MessageInterface::out("CgSolver (Warning): Failure to converge in  ITMAX function evaluations.!!", MessageInterface::DEBUG_ONLY);
      break;
    case 502:
      MessageInterface::out("CgSolver (Warning): Function does not change sign at the endpoints!!", MessageInterface::DEBUG_ONLY);
      break;
    case 503:
      MessageInterface::out("CgSolver (Warning): Red-black indexing is not possible!!", MessageInterface::DEBUG_ONLY);
      break;
    case 601:
      MessageInterface::out("CgSolver (Warning): Successive iterates are not monotone increasing!!", MessageInterface::DEBUG_ONLY);
      break;
    }

    return;
  }
  
}
