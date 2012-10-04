#include <spline/utilities.hpp>

namespace imaging
{
  void interpolation_matrix(const Bspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix)
  {
    spline_utilities_impl::interpolation_matrix(spline, values, matrix);
  }
  
  void interpolation_matrix(const PeriodicBspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix)
  {
    spline_utilities_impl::interpolation_matrix(spline, values, matrix);
  }
  
  void basis_spline_matrix(const PeriodicBspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix) {}
  
  void basis_spline_matrix(const Bspline<float_t> & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix) {}
  
  void basis_spline(size_t i, Bspline<float_t> & spline)
  {
    spline_utilities_impl::basis_spline(i, spline);
  }
  
  void basis_spline(size_t i, PeriodicBspline<float_t> & spline)
  {
    spline_utilities_impl::basis_spline(i, spline);
  }
 
  namespace spline_utilities_impl
  {
    template <class spline_t>
    void basis_spline(size_t i, spline_t & spline)
    {
      for(size_t j = 0; j < spline.n_coefficients(); ++j)
        spline.set_coefficient(j, float_t(0.0));
  
      if(0 <= i && i < spline.n_coefficients())
        spline.set_coefficient(i, float_t(1.0));
    }
    
    template <class spline_t>
    void interpolation_matrix(const spline_t & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix)
    {
      if(spline.n_coefficients() == 0)
        return;
        
      matrix.resize(spline.n_coefficients(), spline.n_coefficients(), false);
      spline_t basis_spline_i(spline), basis_spline_j(spline);
      
      for(size_t i = 0; i < spline.n_coefficients(); ++i)
      {
        basis_spline(i, basis_spline_i);
        for(size_t j = 0; j < spline.n_coefficients(); ++j)
        {
          basis_spline(j, basis_spline_j);
          for(size_t k = 0; k < values.size(); ++k)
            //if(spline.is_in_basis_spline_support(i, values(k)) && spline.is_in_basis_spline_support(j, values(k)))
              matrix(i, j) += basis_spline_i(values(k)) * basis_spline_j(values(k));
        }
      }
    }
    
    template <class spline_t>
    void basis_spline_matrix(const spline_t & spline, const ublas::vector<float_t> & values, ublas::compressed_matrix< float_t > & matrix)
    {
    }
  }
}
