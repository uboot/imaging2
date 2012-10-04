#include <core/matrix_utilities.hpp>

#include <core/utilities.hpp>

namespace imaging
{
  float_t determinant(const ublas::fixed_matrix<float_t, 2, 2> & A)
  {
    return A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
  }

  float_t determinant(const ublas::fixed_matrix<float_t, 3, 3> & A)
  {
    return A(0, 0) * A(1, 1) * A(2, 2) - A(0, 0) * A(1, 2) * A(2, 1) - A(1, 0) * A(0, 1) * A(2, 2)
          +A(1, 0) * A(0, 2) * A(2, 1) + A(2, 0) * A(0, 1) * A(1, 2) - A(2, 0) * A(0, 2) * A(1, 1);
  }


  ublas::fixed_matrix<float_t, 2, 2> inverse(const ublas::fixed_matrix<float_t, 2, 2> & A)
  {
    ublas::fixed_matrix<float_t, 2, 2> C;
    float_t det = determinant(A);
    if(det == 0.0) throw MathException("MathException: Singular matrix in tla::matrix::inverse().");

    C(0, 0) = A(1, 1);
    C(1, 1) = A(0, 0);
    C(0, 1) = - A(0, 1);
    C(1, 0) = - A(1, 0);

    C /= det;

    return C;
  }

  ublas::fixed_matrix<float_t, 3, 3> inverse(const ublas::fixed_matrix<float_t, 3, 3> & A)
  {
    ublas::fixed_matrix<float_t, 3, 3> C;
    float_t det = determinant(A);
    if(det == 0.0) throw MathException("MathException: Singular matrix in tla::matrix::inverse().");
    
    
    C(0, 0) = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
    C(0, 1) = A(0, 2) * A(2, 1) - A(2, 2) * A(0, 1);
    C(0, 2) = A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2);

    C(1, 0) = A(1, 2) * A(2, 0) - A(2, 2) *A(1, 0);
    C(1, 1) = A(0, 0) * A(2, 2) - A(0, 2) *A(2, 0);
    C(1, 2) = A(0, 2) * A(1, 0) - A(1, 2) *A(0, 0);

    C(2, 0) = A(1,0) * A(2,1) - A(2,0) *A(1,1);
    C(2, 1) = A(0,1) * A(2,0) - A(2,1) *A(0,0);
    C(2, 2) = A(0,0) * A(1,1) - A(0,1) *A(1,0);

    C /= det;

    return C;
   }


  ublas::fixed_matrix<float_t, 1, 1> inverse(const ublas::fixed_matrix<float_t, 1, 1> & A)
  {
    ublas::fixed_matrix<float_t, 1, 1> C;
    if(A(0, 0) == 0.0) throw MathException("MathException: Singular matrix in tla::matrix::inverse().");

    C(0, 0) = 1.0 / A(0, 0);

    return C;
  }

  ublas::fixed_matrix<float_t, 2, 2> rotation_matrix(float_t alpha)
  {
    ublas::fixed_matrix<float_t, 2, 2> R;
    
    R(0, 0) = cos(alpha);
    R(0, 1) = -sin(alpha);
    R(1, 0) = -R(0, 1);
    R(1, 1) = R(0, 0);

    return R;
  }
}

