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

#ifndef STATISTIC_LINEARPCA_H
#define STATISTIC_LINEARPCA_H

#include <core/imaging2.hpp>


namespace imaging
{
  /** \ingroup statistic
      \brief Linear principal component analysis.
      
      This class computes a (linear) principal component analysis (PCA) of given data. It provides functions to retrieve the standard deviations of the data within the principal components and to compute the PCA coefficients of vectors and vice-versa.
      
      A LinearPCA object is always set to some current data and a current dimension, either by passing the data and the dimension to the constructor or by calling set_data(). All member functions refer to the current data at the current dimension. The default dimension is the data dimension, i.e. the number of columns of the data matrix.
  */
  class LinearPca
  {
    size_t _dimension;
    size_t _data_dimension;
    
    ublas::vector<float_t> _mean;
    ublas::vector<float_t> _standard_deviations;
    ublas::matrix<float_t> _root_of_covariance;
    ublas::matrix<float_t> _inverse_root_of_covariance;

    void set_data(const ublas::matrix<float_t> & data);
    void set_data(const ublas::matrix<float_t> & data, size_t dimension);
    
  public:
    LinearPca() : _dimension(0), _data_dimension(0) {}
    
    /** Constructs a LinearPca and sets the current data to \em data. The rows of \em data correspond to sample vectors, its columns to the components of the sample vectors. This function results in a PCA which preserves the dimension of the data, i.e. the PCA coefficients of a vector are of the same dimension as the original vector. The current dimension is set to the data dimension.
        The PCA is always computed with respect to the mean of the data. You do not have to center yourself.
    
        Keep in mind that this function computes the PCA, which involves the computation of the data covariance matrix and its eigensystem. */ 
    LinearPca(const ublas::matrix<float_t> & data);
    
    /** Constructs a LinearPca and sets the current data to \em data. The rows of \em data correspond to sample vectors, its columns to the components of the sample vectors. This function results in a PCA which reduces the dimension of the data to \em dimension, i.e. the PCA coefficients of a vector have less components than the original vector. The current dimension is set to \em dimension.
    
    Keep in mind that this function computes the PCA, which involves the computation of the data covariance matrix and its eigensystem. */  
    LinearPca(const ublas::matrix<float_t> & data, size_t dimension);
    
    /** Returns the mean of the current data. */
    const ublas::vector<float_t> & mean() const;
    
    /** Computes a vector from PCA coefficients. The size of \em coefficients must be the same as the current dimension of the PCA. */
    void compute_vector(const ublas::vector<float_t> & coefficients, ublas::vector<float_t> & vector) const;
    
    /** Computes the 2-norm of the PCA coefficients of \em vector. This is the same as the Mahalanobis distance between \em vector and the mean of the current data. */
    float_t norm(const ublas::vector<float_t> & vector) const;
    
    /** Computes the PCA coefficients of \em vector. The size of \em vector must be the same as the dimension of the current data. */
    void compute_coefficients(const ublas::vector<float_t> & vector, ublas::vector<float_t> & coefficients) const;
    
    /** Computes the PCA coefficients of the rows of \em matrix and stores them in the rows of \em coefficients. The number of columns of \em matrix must be the same as the dimension of the current data. */
    void compute_coefficients(const ublas::matrix<float_t> & matrix, ublas::matrix<float_t> & coefficients) const;
    
    /** Returns a reference to the standard deviations within the principal components of the current data. The standard deviations are sorted from the largest to the smallest value. */
    const ublas::vector<float_t> & standard_deviations() const;
    
    /** Returns the dimension of the PCA, i.e. the number of principal components. In case the dimension is smaller than the data dimension, the PCA reduces the dimension of the original data. */
    size_t dimension() const;
    
    /** Returns the dimension of the current data regardless of the actual dimension of the PCA. */ 
    size_t data_dimension() const;
    
    /** Computes the directional derivative into \em direction of the sum of the squared PCA coefficients of \em vector. This is the same as the directional derivative into \em direction of the Mahalanobis distance between \em vector and the mean of the current data . */
    float_t squared_norm_derivative(const ublas::vector<float_t> & vector, const ublas::vector<float_t> & direction) const;
  };
}

#endif
