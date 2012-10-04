#include <shape/ShapeStatistics.hpp>

namespace imaging
{
  namespace shape_statistics_impl
  {
    void constructor(const std::vector<const ShapeInterface*> shape_ptrs, ShapeInterface & mean_shape, ublas::matrix<float_t> & covariance)
    {
      if(shape_ptrs.size() == 0)
        throw Exception("Exception: no shapes in ShapeStatistics::ShapeStatistics()");
      
      const ShapeInterface & reference_shape = *shape_ptrs[0];
      ublas::vector<float_t> current_logarithm(reference_shape.dimension());
      ublas::vector<float_t> mean(reference_shape.dimension());
      
      mean = ublas::scalar_vector<float_t>(mean.size(), 0.0);
      
      for(std::size_t i = 1; i < shape_ptrs.size(); ++i)
      {
        reference_shape.logarithm(*(shape_ptrs[i]), current_logarithm);
        
        mean += current_logarithm;
      }
      
      mean /= float_t(shape_ptrs.size());
       
      reference_shape.exponential(mean, mean_shape);
            
      covariance.resize(reference_shape.dimension(), reference_shape.dimension());
      
      covariance = ublas::scalar_matrix<float_t>(covariance.size1(), covariance.size2(), 0.0);
      
      for(std::size_t i = 0; i < shape_ptrs.size(); ++i)
      {
        mean_shape.logarithm(*(shape_ptrs[i]), current_logarithm);
        
        covariance += outer_prod(current_logarithm, current_logarithm);
      }
      
      // maximum likelihood estimator of the covariance matrix
      covariance /= float_t(shape_ptrs.size());

    }
  
    void compute_statistics(const std::map<std::size_t, float_t> & manual_modes, const ublas::matrix<float_t> & covariance, ublas::matrix<float_t> & root_of_covariance)
    {
      std::size_t n_active_modes = covariance.size1() - manual_modes.size();
    
      ublas::matrix<float_t> temp_covariance(n_active_modes, n_active_modes);
      ublas::matrix<float_t> temp_eigenvectors;
      ublas::vector<float_t> temp_eigenvalues;
      ublas::matrix<float_t> temp_root_of_covariance(n_active_modes, n_active_modes);
      
      for(std::size_t i = 0, i_temp = 0; i_temp < n_active_modes; ++i, ++i_temp)
      {
        while(manual_modes.count(i) > 0) ++i;
        for(std::size_t j = 0, j_temp = 0; j_temp < n_active_modes; ++j, ++j_temp)
        {
          while(manual_modes.count(j) > 0) ++j;
          temp_covariance(i_temp, j_temp) = covariance(i, j);
        }
      }
      
      symmetric_square_root(temp_covariance, temp_root_of_covariance);
      
      root_of_covariance.resize(covariance.size1(), covariance.size2());
      root_of_covariance = ublas::scalar_matrix<float_t>(covariance.size1(), covariance.size2(), 0.0);
      
      for(std::size_t i = 0, i_temp = 0; i_temp < n_active_modes; ++i, ++i_temp)
      {
        while(manual_modes.count(i) > 0) ++i;
        for(std::size_t j = 0, j_temp = 0; j_temp < n_active_modes; ++j, ++j_temp)
        {
          while(manual_modes.count(j) > 0) ++j;
          root_of_covariance(i, j) = temp_root_of_covariance(i_temp, j_temp);
        }
      }
     
      for(std::map<std::size_t, float_t>::const_iterator iter = manual_modes.begin(); iter != manual_modes.end(); ++iter)
        root_of_covariance(iter->first, iter->first) = iter->second;
    }
  }
}

