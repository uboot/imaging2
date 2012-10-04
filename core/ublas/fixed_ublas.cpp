#include <iostream>

typedef double my_float;

#include <core/ublas/fixed_vector.hpp>
#include <core/ublas/fixed_matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/io.hpp>

int main (int argc, char * const argv[]) {

    my_float t = 2.0;
    boost::numeric::ublas::fixed_vector<my_float, 2> u(1.0, 2.0), v;
    v.assign(3.0, 4.0);
    
    std::cout << norm_1(u) << std::endl;
    std::cout << norm_2(u) << std::endl;
    std::cout << norm_inf(u) << std::endl;
    std::cout << - u << std::endl;
    std::cout << u * t << std::endl;
    std::cout << t * u << std::endl;
    std::cout << u + v << std::endl;
    std::cout << u - u << std::endl;
    
    std::cout << inner_prod(u, v) << std::endl;
    std::cout << outer_prod(u, v) << std::endl;
    std::cout << element_prod(u, v) << std::endl;
    std::cout << element_div(u, v) << std::endl;
    
    boost::numeric::ublas::fixed_matrix<my_float, 2, 3> A;
    A = boost::numeric::ublas::scalar_matrix<my_float>(2, 3, 1.0);
    
    boost::numeric::ublas::fixed_matrix<my_float, 2, 3> B;
    B(0, 0) = 1.0;
    B(0, 1) = 2.0;
    B(0, 2) = 3.0;
    B(1, 0) = 4.0;
    B(1, 1) = 5.0;
    B(1, 2) = 6.0;
    
    std::cout << A << std::endl;
    std::cout << B << std::endl;
    
    std::cout << norm_inf(B) << std::endl;
    std::cout << - B << std::endl;
    std::cout << B * t << std::endl;
    std::cout << t * B << std::endl;
    std::cout << A + B << std::endl;
    std::cout << B - B << std::endl;
    
    std::cout << element_prod(A, B) << std::endl;
    std::cout << element_div(A, B) << std::endl;
    
    boost::numeric::ublas::fixed_matrix<my_float, 3, 2> C(trans(B));
    
    std::cout << C << std::endl;
    std::cout << prod(A, C) << std::endl;
    std::cout << prod(trans(C), trans(A)) << std::endl;
    std::cout << prod(C, u) << std::endl;
    std::cout << prod(u, A) << std::endl;
    
    boost::numeric::ublas::banded_matrix<my_float> D(3, 3);
    
    D = 3.0 * boost::numeric::ublas::identity_matrix<my_float>(3);
    
    std::cout << D << std::endl;
    std::cout << prod(D, C) << std::endl;
    std::cout << prod(D, D) << std::endl;
    
    boost::numeric::ublas::matrix_row< boost::numeric::ublas::fixed_matrix<my_float, 2, 3> > row(B, 1);
    boost::numeric::ublas::matrix_column< boost::numeric::ublas::fixed_matrix<my_float, 2, 3> > column(B, 2);
    
    std::cout << row << std::endl;
    std::cout << column << std::endl;
    
    boost::numeric::ublas::fixed_vector<my_float, 3> w(7.0, 7.0, 7.0);
    row = w;
    std::cout << B << std::endl;
    
    column = u;
    std::cout << B << std::endl;
    
}
