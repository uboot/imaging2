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

#ifndef CORE_UBLAS_FIXEDSIZEARRAY_H
#define CORE_UBLAS_FIXEDSIZEARRAY_H

#include <boost/numeric/ublas/storage.hpp>

namespace boost { namespace numeric { namespace ublas {

    // Fixed size array 
    template<class T, std::size_t N> 
    class fixed_size_array: 
        public storage_array<fixed_size_array<T, N> > { 

        typedef fixed_size_array<T, N> self_type; 
    public: 
        typedef std::size_t size_type; 
        typedef std::ptrdiff_t difference_type; 
        typedef T value_type; 
        typedef const T &const_reference; 
        typedef T &reference; 
        typedef const T *const_pointer; 
        typedef T *pointer; 
        typedef const_pointer const_iterator; 
        typedef pointer iterator; 

        // Construction and destruction 
        BOOST_UBLAS_INLINE 
        fixed_size_array () { 
        } 
        explicit BOOST_UBLAS_INLINE 
        fixed_size_array (size_type size) { 
            BOOST_UBLAS_CHECK (size == N, bad_size ()); 
        }
        BOOST_UBLAS_INLINE 
        fixed_size_array (size_type size, const value_type &init) { 
            BOOST_UBLAS_CHECK (size == N, bad_size ()); 
            std::fill (data_, data_ + N, init) ; 
        } 
        BOOST_UBLAS_INLINE 
        fixed_size_array (const fixed_size_array &c) : storage_array<fixed_size_array<T, N> >() { 
            std::copy (c.data_, c.data_ + N, data_); 
        } 

        // Resizing 
        BOOST_UBLAS_INLINE 
        void resize (size_type size) { 
            BOOST_UBLAS_CHECK (size == N, bad_size ()); 
        } 
        BOOST_UBLAS_INLINE 
        void resize (size_type size, value_type init) { 
            BOOST_UBLAS_CHECK (size == N, bad_size ()); 
            std::fill (data_, data_ + N, init); 
        }

        // Random Access Container 
        BOOST_UBLAS_INLINE 
        size_type max_size () const { 
            return N; 
        } 

        BOOST_UBLAS_INLINE 
        bool empty () const { 
            return N == 0; 
        } 

        BOOST_UBLAS_INLINE 
        size_type size () const { 
            return N; 
        } 

        // Element access 
        BOOST_UBLAS_INLINE 
        const_reference operator [] (size_type i) const { 
            BOOST_UBLAS_CHECK (i < N, bad_index ()); 
            return data_ [i]; 
        } 
        BOOST_UBLAS_INLINE 
        reference operator [] (size_type i) { 
            BOOST_UBLAS_CHECK (i < N, bad_index ()); 
            return data_ [i]; 
        } 

        // Assignment 
        BOOST_UBLAS_INLINE 
        fixed_size_array &operator = (const fixed_size_array &a) { 
            if (this != &a) { 
                std::copy (a.data_, a.data_ + N, data_); 
            } 
            return *this; 
        } 
        BOOST_UBLAS_INLINE 
        fixed_size_array &assign_temporary (fixed_size_array &a) { 
            *this = a; 
            return *this; 
        } 

        // Swapping 
        BOOST_UBLAS_INLINE 
        void swap (fixed_size_array &a) { 
            if (this != &a) { 
                std::swap_ranges (data_, data_ + N, a.data_); 
            } 
        } 
        BOOST_UBLAS_INLINE 
        friend void swap (fixed_size_array &a1, fixed_size_array &a2) { 
            a1.swap (a2); 
        } 

        BOOST_UBLAS_INLINE 
        const_iterator begin () const { 
            return data_; 
        } 
        BOOST_UBLAS_INLINE 
        const_iterator end () const { 
            return data_ + N; 
        } 

        BOOST_UBLAS_INLINE 
        iterator begin () { 
            return data_; 
        } 
        BOOST_UBLAS_INLINE 
        iterator end () { 
            return data_ + N; 
        } 

        // Reverse iterators 
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator; 
        typedef std::reverse_iterator<iterator> reverse_iterator; 

        BOOST_UBLAS_INLINE 
        const_reverse_iterator rbegin () const { 
            return const_reverse_iterator (end ()); 
        } 
        BOOST_UBLAS_INLINE 
        const_reverse_iterator rend () const { 
            return const_reverse_iterator (begin ()); 
        } 
        BOOST_UBLAS_INLINE 
        reverse_iterator rbegin () { 
            return reverse_iterator (end ()); 
        } 
        BOOST_UBLAS_INLINE 
        reverse_iterator rend () { 
            return reverse_iterator (begin ()); 
        } 

    private: 
        value_type data_ [N]; 
    };

    // Fixed size array of lenth zero
    template<class T> 
    class fixed_size_array<T, 0> : 
        public storage_array<fixed_size_array<T, 0> > { 

        typedef fixed_size_array<T, 0> self_type; 
    public: 
        typedef std::size_t size_type; 
        typedef std::ptrdiff_t difference_type; 
        typedef T value_type; 
        typedef const T &const_reference; 
        typedef T &reference; 
        typedef const T *const_pointer; 
        typedef T *pointer; 
        typedef const_pointer const_iterator; 
        typedef pointer iterator; 

        // Construction and destruction 
        BOOST_UBLAS_INLINE 
        fixed_size_array () { 
        } 
        explicit BOOST_UBLAS_INLINE 
        fixed_size_array (size_type size) { 
            BOOST_UBLAS_CHECK (size == 0, bad_size ()); 
        } 
        BOOST_UBLAS_INLINE 
        fixed_size_array (size_type size, const value_type &init) { 
            BOOST_UBLAS_CHECK (size == 0, bad_size ()); 
        } 
        BOOST_UBLAS_INLINE 
        fixed_size_array (const fixed_size_array &c) { 
        } 

        // Resizing 
        BOOST_UBLAS_INLINE 
        void resize (size_type size) { 
            BOOST_UBLAS_CHECK (size == 0, bad_size ()); 
        } 
        BOOST_UBLAS_INLINE 
        void resize (size_type size, value_type init) { 
            BOOST_UBLAS_CHECK (size == 0, bad_size ()); 
        } 

        // Random Access Container 
        BOOST_UBLAS_INLINE 
        size_type max_size () const { 
            return 0; 
        } 

        BOOST_UBLAS_INLINE 
        bool empty () const { 
            return 0 == 0; 
        } 

        BOOST_UBLAS_INLINE 
        size_type size () const { 
            return 0; 
        } 

        // Assignment 
        BOOST_UBLAS_INLINE 
        fixed_size_array &operator = (const fixed_size_array &a) { 
            return *this; 
        } 
        BOOST_UBLAS_INLINE 
        fixed_size_array &assign_temporary (fixed_size_array &a) { 
            return *this; 
        } 

        // Swapping 
        BOOST_UBLAS_INLINE 
        void swap (fixed_size_array &a) { 
        } 
        BOOST_UBLAS_INLINE 
        friend void swap (fixed_size_array &a1, fixed_size_array &a2) { 
        } 

        BOOST_UBLAS_INLINE 
        const_iterator begin () const { 
            return 0; 
        } 
        BOOST_UBLAS_INLINE 
        const_iterator end () const { 
            return 0; 
        } 

        BOOST_UBLAS_INLINE 
        iterator begin () { 
            return 0; 
        } 
        BOOST_UBLAS_INLINE 
        iterator end () { 
            return 0; 
        } 

        // Reverse iterators 
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator; 
        typedef std::reverse_iterator<iterator> reverse_iterator; 

        BOOST_UBLAS_INLINE 
        const_reverse_iterator rbegin () const { 
            return const_reverse_iterator (end ()); 
        } 
        BOOST_UBLAS_INLINE 
        const_reverse_iterator rend () const { 
            return const_reverse_iterator (begin ()); 
        } 
        BOOST_UBLAS_INLINE 
        reverse_iterator rbegin () { 
            return reverse_iterator (end ()); 
        } 
        BOOST_UBLAS_INLINE 
        reverse_iterator rend () { 
            return reverse_iterator (begin ()); 
        } 
    }; 
}}}

#endif
