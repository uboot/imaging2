/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CORE_UBLAS_FIXEDVECTOR_H
#define CORE_UBLAS_FIXEDVECTOR_H

#include <boost/numeric/ublas/vector.hpp>
#include <core/ublas/fixed_size_array.hpp>

namespace boost { namespace numeric { namespace ublas {

    /** \ingroup core
        \brief uBLAS compatible class template for vectors of fixed size.
    
        This class implements a uBLAS compatible vector class. The template parameters \em N specifies the size of the vector which has to be known at compile time. See <a href="http://www.boost.org/libs/numeric/ublas/doc/index.htm">uBLAS documentation</a> for details.
    */
    template<class T, std::size_t N>
    class fixed_vector:
        public vector<T, fixed_size_array<T, N> > {

        typedef vector<T, fixed_size_array<T, N> > vector_type;
    public:
    /** \cond */
        typedef typename vector_type::size_type size_type;
        static const size_type max_size = N;
        static const size_type dimension = N;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        fixed_vector ():
            vector_type (N) {}
        /*BOOST_UBLAS_INLINE
        fixed_vector (size_type size):
            vector_type (N) {}*/
        BOOST_UBLAS_INLINE
        fixed_vector (const fixed_vector &v):
            vector_type (v) {}
        template<class A2>              // Allow vector<T,bounded_array<N> construction
        BOOST_UBLAS_INLINE
        fixed_vector (const vector<T, A2> &v):
            vector_type (v) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        fixed_vector (const vector_expression<AE> &ae):
            vector_type (ae) {}
        BOOST_UBLAS_INLINE
        ~fixed_vector () {}

        // Assignment
        BOOST_UBLAS_INLINE
        fixed_vector &operator = (const fixed_vector &v) {
            vector_type::operator = (v);
            return *this;
        }
        template<class A2>         // Generic vector assignment
        BOOST_UBLAS_INLINE
        fixed_vector &operator = (const vector<T, A2> &v) {
            vector_type::operator = (v);
            return *this;
        }
        template<class C>          // Container assignment without temporary
        BOOST_UBLAS_INLINE
        fixed_vector &operator = (const vector_container<C> &v) {
            vector_type::operator = (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        fixed_vector &operator = (const vector_expression<AE> &ae) {
            vector_type::operator = (ae);
            return *this;
        }
        
        
        // componentwise initialization
    
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v)
        { 
          assign(v);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1) 
        {
            assign(v0, v1);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2) 
        {
            assign(v0, v1, v2);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2, const T & v3) 
        {
            assign(v0, v1, v2, v3);
        }
          
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4) 
        {
          assign(v0, v1, v2, v3, v4);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5) 
        {
          assign(v0, v1, v2, v3, v4, v5);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6) 
        {
          assign(v0, v1, v2, v3, v4, v5, v6);
        }
        
        BOOST_UBLAS_INLINE
        explicit fixed_vector(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6, const T & v7) 
        {
          assign(v0, v1, v2, v3, v4, v5, v6, v7);
        }
        
        template<class DATA_t>
        BOOST_UBLAS_INLINE
        explicit fixed_vector (const fixed_vector<DATA_t, N> &v)
        {
            for(typename vector<T, fixed_size_array<T, N> >::iterator iter = vector<T, fixed_size_array<T, N> >::begin(); iter != vector<T, fixed_size_array<T, N> >::end(); ++iter)
                  *iter = T(v(iter.index()));
        }
        
        BOOST_UBLAS_INLINE
        void assign(const T & v) 
        {
            for(typename vector<T, fixed_size_array<T, N> >::iterator iter = vector<T, fixed_size_array<T, N> >::begin(); iter != vector<T, fixed_size_array<T, N> >::end(); ++iter)
                  *iter = v;
        }

        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1)
        {
            (*this)(0) = v0;
            (*this)(1) = v1;
        }

        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2) 
        {
            (*this)(0) = v0;
            (*this)(1) = v1;
            (*this)(2) = v2;
        }

        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2, const T & v3) 
        {
            (*this)(0) = v0;
            (*this)(1) = v1;
            (*this)(2) = v2;
            (*this)(3) = v3;
        }
          
        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4) 
        {
          (*this)(0) = v0;
          (*this)(1) = v1;
          (*this)(2) = v2;
          (*this)(3) = v3;
          (*this)(4) = v4;
        }
        
        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5) 
        {
          (*this)(0) = v0;
          (*this)(1) = v1;
          (*this)(2) = v2;
          (*this)(3) = v3;
          (*this)(4) = v4;
          (*this)(5) = v5;
        }
        
        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6) 
        {
          (*this)(0) = v0;
          (*this)(1) = v1;
          (*this)(2) = v2;
          (*this)(3) = v3;
          (*this)(4) = v4;
          (*this)(5) = v5;
          (*this)(6) = v6;
        }
        
        BOOST_UBLAS_INLINE
        void assign(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6, const T & v7) 
        {
          (*this)(0) = v0;
          (*this)(1) = v1;
          (*this)(2) = v2;
          (*this)(3) = v3;
          (*this)(4) = v4;
          (*this)(5) = v5;
          (*this)(6) = v6;
          (*this)(7) = v7;
        }
    /** \endcond */
    };
    
    // comparision
    template<class T, std::size_t N>
    bool operator==(const fixed_vector<T, N> & lhs, const fixed_vector<T, N> & rhs)
    {
          for(typename fixed_vector<T, N>::const_iterator iter = lhs.begin(); iter != lhs.end(); ++iter)
                if(*iter != rhs(iter.index())) return false;
          
          return true;
    }
    
    // comparision
    template<class T, std::size_t N>
    bool operator!=(const fixed_vector<T, N> & lhs, const fixed_vector<T, N> & rhs)
    {
          return ! (lhs == rhs);
    }
}}}

#endif
