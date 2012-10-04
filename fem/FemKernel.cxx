#include <fem/FemKernel.hpp>


namespace imaging
{
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 1, 0> & derivative)
  {
    return 1.0;
  }

  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 1, 1> & derivative)
  {
    return fabs(derivative(0, 0));
  }
  
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 2, 1> & derivative)
  {
    return sqrt(derivative(0,0) * derivative(0,0) + derivative(1,0) * derivative(1,0));
  }
  
  template <>
  float_t transform_det(const ublas::fixed_matrix<float_t, 3, 2> & A)
  {
    return sqrt(A(0, 0)*A(0,0) * A(1, 1)* A(1, 1)+ 
                A(0, 0)*A(0, 0)*A(2, 1)*A(2, 1)+
                A(1, 0)*A(1, 0)*A(0, 1)*A(0, 1)+
		A(1, 0)*A(1, 0)*A(2, 1)*A(2, 1)+
	        A(2, 0)*A(2, 0)*A(0, 1)*A(0, 1)+
                A(2, 0)*A(2, 0)*A(1, 1)*A(1, 1)
		-1*A(0, 0)*A(0, 1)*A(1, 0)*A(1, 1)
		-1*A(0, 0)*A(0, 1)*A(2, 0)*A(2, 1)
		-1*A(1, 0)*A(1, 1)*A(2, 0)*A(2, 1));
  }

}


