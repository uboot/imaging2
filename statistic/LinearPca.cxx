#include <statistic/LinearPca.hpp>

#include <lapack/linear_algebra.hpp>

namespace imaging
{
  LinearPca::LinearPca(const ublas::matrix<float_t> & data)
  {
    set_data(data, data.size2());
  }
  
  LinearPca::LinearPca(const ublas::matrix<float_t> & data, size_t dimension)
  {
    if(dimension > data.size2())
      throw Exception("Specified too large dimension in LinearPca::LinearPca().");
    
    set_data(data, dimension);
  }
  
  void LinearPca::set_data(const ublas::matrix<float_t> & data)
  {
    set_data(data, data.size2());
  }
  
  void LinearPca::set_data(const ublas::matrix<float_t> & data, size_t dimension)
  {
    if(dimension > data.size2())
      throw Exception("Specified too large dimension in LinearPca::set_data().");
    
    _dimension = dimension;
    _data_dimension = data.size2();
    
    size_t n_samples = data.size1();

    _mean.resize(_data_dimension);
    _mean = ublas::scalar_vector<float_t>(_data_dimension, 0.0);
        
    for(std::size_t i = 0; i < n_samples; ++i)
      _mean += ublas::matrix_row< const ublas::matrix<float_t> >(data, i);
    
    _mean /= float_t(n_samples);
      
    ublas::matrix<float_t> covariance(_data_dimension, _data_dimension);
    
    covariance = ublas::scalar_matrix<float_t>(_data_dimension, _data_dimension, 0.0);
    
    for(std::size_t i = 0; i < n_samples; ++i)
      covariance += outer_prod(ublas::matrix_row< const ublas::matrix<float_t> >(data, i) - _mean,
                               ublas::matrix_row< const ublas::matrix<float_t> >(data, i) - _mean);
    
    
    // maximum likelihood estimator of the covariance matrix
    covariance /= float_t(n_samples);

    _root_of_covariance.resize(_data_dimension, _dimension);
    _inverse_root_of_covariance.resize(_dimension, _data_dimension);
    _standard_deviations.resize(_dimension);
    
    ublas::matrix<float_t> eigenvectors(_data_dimension, _data_dimension);
    ublas::vector<float_t> eigenvalues(_data_dimension);
    
    eigensystem(covariance, eigenvectors, eigenvalues);
    
    for(size_t i = 0; i < _dimension; ++i)
      _standard_deviations(i) = sqrt(fabs(eigenvalues(i)));

    for(ublas::matrix<float_t>::iterator1 iter1 = _root_of_covariance.begin1(); iter1 != _root_of_covariance.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = eigenvectors(iter2.index1(), iter2.index2()) * 
                 _standard_deviations(iter2.index2());
    
    for(ublas::matrix<float_t>::iterator1 iter1 = _inverse_root_of_covariance.begin1(); iter1 != _inverse_root_of_covariance.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = eigenvectors(iter2.index2(), iter2.index1()) / 
                 _standard_deviations(iter2.index1());
  }
    
  const ublas::vector<float_t> & LinearPca::mean() const
  {
    return _mean;
  }
  
  void LinearPca::compute_vector(const ublas::vector<float_t> & coefficients, ublas::vector<float_t> & vector) const
  {
    vector.resize(_data_dimension);
    
    vector = prod(_root_of_covariance, coefficients) + _mean;
  }
  
  float_t LinearPca::norm(const ublas::vector<float_t> & vector) const
  {
    ublas::vector<float_t> coefficients;
    compute_coefficients(vector, coefficients);
    return norm_2(coefficients);
  }
  
  void LinearPca::compute_coefficients(const ublas::vector<float_t> & vector, ublas::vector<float_t> & coefficients) const
  {
    coefficients = prod(_inverse_root_of_covariance, vector - _mean);
  }
  
  void LinearPca::compute_coefficients(const ublas::matrix<float_t> & matrix, ublas::matrix<float_t> & coefficients) const
  {
    coefficients = trans(prod(_inverse_root_of_covariance, trans(matrix) - outer_prod(_mean, ublas::scalar_vector<float_t>(matrix.size1(), 1.0))));
  }
  
  const ublas::vector<float_t> & LinearPca::standard_deviations() const
  {
    return _standard_deviations;
  }
  
  size_t LinearPca::dimension() const
  {
    return _dimension;
  }
  
  size_t LinearPca::data_dimension() const
  {
    return _data_dimension;
  }
  
  float_t LinearPca::squared_norm_derivative(const ublas::vector<float_t> & vector, const ublas::vector<float_t> & direction) const
  {
    return 2.0 * inner_prod(prod(_inverse_root_of_covariance, vector - _mean), prod(_inverse_root_of_covariance, direction));
  }
}
