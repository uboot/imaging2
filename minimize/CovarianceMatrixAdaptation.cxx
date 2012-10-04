#include <minimize/CovarianceMatrixAdaptation.hpp>

#include <core/utilities.hpp>
#include <core/MessageInterface.hpp>
#include <core/distribution_utilities.hpp>
#include <lapack/linear_algebra.hpp>

#include <map>

namespace imaging
{
  CovarianceMatrixAdaptation::CovarianceMatrixAdaptation(EnergyInterface & energy, float_t sigma, float_t min_update) :
      _energy(energy)
  {
    size_t lambda = 4 + size_t(floor(3.0 * log(float_t(_energy.dimension()))));
    init(sigma, min_update, lambda);
  }
  
  
  CovarianceMatrixAdaptation::CovarianceMatrixAdaptation(EnergyInterface & energy, float_t sigma, float_t min_update, std::size_t lambda) :
      _energy(energy)
  {
    init(sigma, min_update, lambda);
  }
      
  void CovarianceMatrixAdaptation::init(float_t sigma, float_t min_update, std::size_t lambda)
  { 
    _lambda = lambda;
    _sigma = sigma;
    _mu = _lambda / 2;
    _w.resize(_mu);
    _n = _energy.dimension();
    _C.resize(_n, _n);
    _p_sigma.resize(_n);
    _p_c.resize(_n);
    _min_update = min_update;
    _terminated = false;
      
    ublas::vector<float_t> w_prime(_mu);
    
    float_t mu_prime = (float_t(_lambda) - 1.0) / 2.0;
    
    for(std::size_t i = 0; i < _mu; ++i)
      w_prime(i) = log(mu_prime + 1.0) - log(float_t(i) + 1.0);
      
    float_t sum_of_w_prime = inner_prod(w_prime, ublas::scalar_vector<float_t>(w_prime.size(), 1.0));
      
    _w = w_prime / sum_of_w_prime;
      
    _mu_eff = 1.0 / inner_prod(_w, _w);
      
    _c_sigma = (_mu_eff + 2.0) / (float_t(_n) + _mu_eff + 3.0);
    
    _d_sigma = 1.0 + 2.0 * max(0.0, sqrt( (_mu_eff - 1.0) / (float_t(_n) + 1.0) ) - 1.0 ) + _c_sigma;
    
    _c_c = 4.0 / (float_t(_n) + 4.0);
    
    _mu_cov = _mu_eff;
    
    _c_cov = 1.0 / _mu_cov * 2.0 / square(float_t(_n) + SQUARE_ROOT_2) +
            (1.0 - 1.0 / _mu_cov) * min(1.0, (2 * _mu_eff - 1) / (square(float_t(_n) + 2) + _mu_eff) );
            
    
    _C = ublas::identity_matrix<float_t>(_n);
    
    _p_sigma = ublas::scalar_vector<float_t>(_n, 0.0);
    
    _p_c = ublas::scalar_vector<float_t>(_n, 0.0);
  }
    

  bool CovarianceMatrixAdaptation::minimize(std::size_t n_max_steps, std::size_t & n_steps)
  {
    if(_terminated)
    {
      n_steps = 0;
      return true;
    }
    
    float_t expected_norm = SQUARE_ROOT_2 * tgamma(float_t(_n + 1.0) / 2.0) / tgamma(float_t(_n) / 2.0); 
    ublas::vector<float_t> m(_n); 
    ublas::vector<float_t> y_w(_n);
    ublas::matrix<float_t> sum_outer_product_y;
    ublas::matrix<float_t> BD(_n, _n);
    ublas::matrix<float_t> inverse_square_root_C(_n, _n);
    
    std::vector< ublas::vector<float_t> > z(_lambda);
    std::vector< ublas::vector<float_t> > y(_lambda);
    std::vector< ublas::vector<float_t> > x(_lambda);
    
    for(std::size_t i = 0; i < _lambda; ++i) z[i].resize(_n);
    for(std::size_t i = 0; i < _lambda; ++i) y[i].resize(_n);
    for(std::size_t i = 0; i < _lambda; ++i) x[i].resize(_n);
      
    std::multimap<float_t, std::size_t> energies;
    
    m = _energy.current_argument();
    
    for(n_steps = 0; n_steps < n_max_steps; ++n_steps)
    {      
      MessageInterface::out("CovarianceMatrixAdaptation: Step " + boost::lexical_cast<std::string>(n_steps + 1), MessageInterface::DEBUG_ONLY);
      
      symmetric_square_root(_C, BD);
      
      for(std::size_t k = 0; k < _lambda; ++k)
        for(std::size_t i = 0; i < _n; ++i)
          z[k](i) = normal_distribution();
          
      for(std::size_t k = 0; k < _lambda; ++k)
        y[k] = prod(BD, z[k]);
          
      for(std::size_t k = 0; k < _lambda; ++k)
        x[k] = m + _sigma * y[k];
        
      energies.clear();
      
      for(std::size_t k = 0; k < _lambda; ++k)
      {
        _energy.current_argument() = x[k];
          
        // the actual computational effort is hidden here:
        _energy.set_argument();
        
        energies.insert(std::pair<float_t, std::size_t>(_energy.current_energy(), k));
      }
      
      MessageInterface::out("Current minimal functional value: " + boost::lexical_cast<std::string>((*energies.begin()).first), MessageInterface::DEBUG_ONLY, +1);
      
      y_w = ublas::scalar_vector<float_t>(_n, 0.0);
      
      std::size_t i = 0;
      for(std::multimap<float_t, std::size_t>::const_iterator iter = energies.begin();
          iter != energies.end() && i < _mu; ++iter, ++i)
        y_w += _w(i) * y[iter->second];
          
      m += _sigma * y_w;
      
      inverse_square_root(_C, inverse_square_root_C);
      
      _p_sigma = (1 - _c_sigma) * _p_sigma + sqrt(_c_sigma * (2 - _c_sigma) * _mu_eff) * prod(inverse_square_root_C, y_w);
  
      _sigma = exp( _c_sigma / _d_sigma * ( norm_2(_p_sigma) / expected_norm - 1 ) );
      
      float_t h_sigma;
      if( norm_2(_p_sigma) / sqrt( 1 - pow(1 - _c_sigma, 2.0 * (float_t(n_steps) + 1.0) ) ) <
          ( 1.4 + 2.0 / (float_t(_n) + 1.0) ) * expected_norm )
        h_sigma = 1.0;
      else
        h_sigma = 0.0;
        
      float_t delta_h_sigma = (1 - h_sigma) * _c_c * (2 - _c_c);
      
      _p_c = (1 - _c_c) * _p_c + h_sigma * sqrt( _c_c * (2.0 - _c_c) * _mu_eff ) * y_w;
      
      sum_outer_product_y = ublas::scalar_matrix<float_t>(_n, _n, 0.0);
      
      i = 0;
      for(std::multimap<float_t, std::size_t>::const_iterator iter = energies.begin();
          iter != energies.end() && i < _mu; ++iter, ++i)
        sum_outer_product_y += _w(i) * outer_prod(y[iter->second], y[iter->second]);
        
      _C = (1.0 - _c_cov) * _C + _c_cov/_mu_cov * ( outer_prod(_p_c, _p_c) + delta_h_sigma * _C ) +
          _c_cov * ( 1.0 - 1.0/_mu_cov ) * sum_outer_product_y;
      
      if(_sigma * norm_frobenius(BD) < _min_update) 
      {
        _terminated = true;
        break;
      }
    }
    
    _energy.current_argument() = m;
    _energy.set_argument();
    
    return _terminated;
  }

}

