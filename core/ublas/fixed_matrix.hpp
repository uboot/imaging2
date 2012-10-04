#ifndef CORE_UBLAS_FIXED_MATRIX_H
#define CORE_UBLAS_FIXED_MATRIX_H

#include <core/ublas/fixed_size_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace boost { namespace numeric { namespace ublas {

    /** \ingroup core
        \brief uBLAS compatible class template for matrices of fixed size.
    
        This class implements a uBLAS compatible matrix class. The template parameters \em M and \em N specify the size of the matrix which have to be known at compile time. See <a href="http://www.boost.org/libs/numeric/ublas/doc/index.htm">uBLAS documentation</a> for details.
    */
    template<class T, std::size_t M, std::size_t N, class L = row_major>
    class fixed_matrix:
        public matrix<T, L, bounded_array<T, M * N> > {

        typedef matrix<T, L, bounded_array<T, M * N> > matrix_type;
    public:
    /** \cond */
        typedef typename matrix_type::size_type size_type;
        static const size_type max_size1 = M;
        static const size_type max_size2 = N;
        static const size_type dimension1 = M;
        static const size_type dimension2 = N;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        fixed_matrix ():
            matrix_type (M, N) {}
        /*BOOST_UBLAS_INLINE
        fixed_matrix (size_type size1, size_type size2):
            matrix_type (size1, size2) {}*/
        BOOST_UBLAS_INLINE
        fixed_matrix (const fixed_matrix &m):
            matrix_type (m) {}
        template<class A2>              // Allow matrix<T, L, fixed_size_array<M,N> > construction
        BOOST_UBLAS_INLINE
        fixed_matrix (const matrix<T, L, A2> &m):
            matrix_type (m) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        fixed_matrix (const matrix_expression<AE> &ae):
            matrix_type (ae) {}
        BOOST_UBLAS_INLINE
        ~fixed_matrix () {}

        // Assignment
        BOOST_UBLAS_INLINE
        fixed_matrix &operator = (const fixed_matrix &m) {
            matrix_type::operator = (m);
            return *this;
        }
        template<class L2, class A2>        // Generic matrix assignment
        BOOST_UBLAS_INLINE
        fixed_matrix &operator = (const matrix<T, L2, A2> &m) {
            matrix_type::operator = (m);
            return *this;
        }
        template<class C>          // Container assignment without temporary
        BOOST_UBLAS_INLINE
        fixed_matrix &operator = (const matrix_container<C> &m) {
            matrix_type::operator = (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        fixed_matrix &operator = (const matrix_expression<AE> &ae) {
            matrix_type::operator = (ae);
            return *this;
        }
        
        
        // Componentwise initialization
        void assign(const T & a) 
        {
            for(typename matrix<T, L, bounded_array<T, M * N> >::iterator1 iter1 = matrix<T, L, bounded_array<T, M * N> >::begin1(); iter1 != matrix<T, L, bounded_array<T, M * N> >::end1(); ++iter1)
            for(typename matrix<T, L, bounded_array<T, M * N> >::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
              *iter2 = a;
        }
    /** \endcond */
    };

}}}

#endif
