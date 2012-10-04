// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_SHAPESTATISTICS_H
#define SHAPE_SHAPESTATISTICS_H

#include <shape/ShapeInterface.hpp>
#include <lapack/linear_algebra.hpp>
#include <core/distribution_utilities.hpp>

#include <map>
#include <boost/shared_ptr.hpp>


namespace imaging
{
  namespace shape_statistics_impl
  {
    void constructor(const std::vector<const ShapeInterface*> shape_ptrs, ShapeInterface & mean_shape, ublas::matrix<float_t> & covariance);
    void compute_statistics(const std::map<std::size_t, float_t> & manual_modes, const ublas::matrix<float_t> & covariance, ublas::matrix<float_t> & root_of_covariance);
  }

  /** \ingroup shape
      \brief Computes the mean shape and the covariance of a shape distribution and samples from this distribution.
      
      This class computes the mean shape and the covariance in the tangent space at the mean shape of a set of shapes. In addition, it is possible to manually set parameters of the distribution (such as position and rotation) to remove dependence of the statistics on them. After the computation of the statistics the user can sample shapes from a normal distribution on the tangent space of the mean shape. The covariance matrix of this distribution is the same as the one of the sample data (with exception of the manually set modes).
      
      The shapes must be elements of a common shape manifold, i.e. \em shape_t must be derived from ShapeInterface and the Riemannian exponential and logarithm of the shapes must be compatible to each other. More precisely, all shapes passed to a ShapeStatistics object must accept all the others as arguments in their exponential and logarithm member functions.
      
      Note that the user must call compute_statistics() before querying a ShapeStatistic object, i.e. before calling mean_shape() or any of the sample members.
      */
  template <class shape_t>
  class ShapeStatistics
  {
    const static std::size_t SHAPE_DIMENSION = shape_t::SHAPE_DIMENSION;
    
    ublas::matrix<float_t> _covariance;
    ublas::matrix<float_t> _root_of_covariance;
    
    shape_t _mean_shape;
    
    std::map<std::size_t, float_t> _manual_modes;
     
  public:
  
    /** Construct a ShapeStatistics object by manually setting its mean. Note that you to call set_mode_manually() of all (dimension of the mean shape) modes if you choose to construct the statistics from the mean shape rather than shape samples. The member compute_statistics() must be called before querying the statistics. */
    ShapeStatistics(const shape_t & mean_shape) : 
      _mean_shape(mean_shape),
      _covariance(ublas::identity_matrix<float_t>(mean_shape.dimension())),
      _root_of_covariance(ublas::identity_matrix<float_t>(mean_shape.dimension())) {}
    
    /** Construct a ShapeStatistics from a set of shapes. The shapes must be elements of a common shape manifold, i.e. \em shape_t must be derived from ShapeInterface and the Riemannian exponential and logarithm of the shapes must be compatible to each other. More precisely, all shapes passed to a ShapeStatistics object must accept all the others as arguments in their exponential and logarithm member functions. In particular, this implies that the dimension of each of the shapes must be the same. The member compute_statistics() must be called before querying the statistics. */
    ShapeStatistics(const std::vector<shape_t> & shapes)
    {
      std::vector<const ShapeInterface *> shape_ptrs(shapes.size());
      for(size_t i = 0; i < shapes.size(); ++i)
        shape_ptrs[i] = &(shapes[i]);
      
      shape_statistics_impl::constructor(shape_ptrs, _mean_shape, _covariance);
    }
    
    /** Construct a ShapeStatistics from a set of shapes. The shapes must be elements of a common shape manifold, i.e. \em shape_t must be derived from ShapeInterface and the Riemannian exponential and logarithm of the shapes must be compatible to each other. More precisely, all shapes passed to a ShapeStatistics object must accept all the others as arguments in their exponential and logarithm member functions. In particular, this implies that the dimension of each of the shapes must be the same. The member compute_statistics() must be called before querying the statistics. */
    ShapeStatistics(const std::vector< boost::shared_ptr<shape_t> > & shapes)
    {
      std::vector<const ShapeInterface *> shape_ptrs(shapes.size());
      for(size_t i = 0; i < shapes.size(); ++i)
        shape_ptrs[i] = shapes[i].get();
      
      shape_statistics_impl::constructor(shape_ptrs, _mean_shape, _covariance);
    }
    
    /** Forces the <em>mode</em>-th diagonal entry of the covariance matrix to be set to \em deviation and sets all other entries in the corresponding row and column of the covariance matrix to zero. The member compute_statistics() must be called before querying the statistics. */
    void set_mode_manually(std::size_t mode, float_t deviation)
    {
      _manual_modes.insert(std::pair<std::size_t, float_t>(mode, deviation));
    }
    
    /** Computes the statistics and takes into account the manually set modes. After calling this function you can query the statistics. */
    void compute_statistics()
    {
      shape_statistics_impl::compute_statistics(_manual_modes, _covariance, _root_of_covariance);
    }
    
    /** Maps the \em coefficients in the PCA space of the sample data to the corresponding tangent \em vector in the tangent space at the mean shape. */
    void shape_vector(const ublas::vector<float_t> & coefficients, ublas::vector<float_t> & vector) const
    {
      vector.resize(dimension());
      vector = prod(_root_of_covariance, coefficients);
    }
    
    /** Maps the \em coefficients in the PCA space of the sample data to the corresponding tangent \em vector in the tangent space at the mean shape and stores the squared Mahalanobis distance of the vector to the mean in \em squared_distance. */
    void shape_vector(const ublas::vector<float_t> & coefficients, ublas::vector<float_t> & vector, float_t & squared_distance) const
    {
      vector.resize(dimension());
      vector = prod(_root_of_covariance, coefficients);
      
      squared_distance = 0.0;
      for(std::size_t i = 0; i < coefficients.size(); ++i)
        if(_manual_modes.count(i) == 0)
          squared_distance += square(coefficients(i));
    }
  
    /** Maps the \em coefficients in the PCA space of the sample data to the corresponding \em shape. */
    void shape_sample(const ublas::vector<float_t> & coefficients, shape_t & out_shape) const
    {
      float_t squared_distance;
      shape_sample(coefficients, out_shape, squared_distance);
    }

    /** Maps the \em coefficients in the PCA space of the sample data to the corresponding \em shape and stores the squared Mahalanobis distance of the shape sample to the mean shape in \em squared_distance. The manually set modes are not reflected in the distance! */
    void shape_sample(const ublas::vector<float_t> & coefficients, shape_t & out_shape, float_t & squared_distance) const
    { 
      ublas::vector<float_t> sample_vector(_mean_shape.dimension());
      shape_vector(coefficients, sample_vector, squared_distance);
  
      _mean_shape.exponential(sample_vector, out_shape);
    }  
    
    /** Maps the \em coefficients in the PCA space of the sample data to gradient of the squared Mahalanobis distance to the mean shape. The manually set modes are not reflected in the distance! */
    void gradient(const ublas::vector<float_t> & coefficients, ublas::vector<float_t> & gradient) const
    { 
      gradient.resize(_mean_shape.dimension());
      
      gradient = 2 * coefficients;
      
      for(std::map<std::size_t, float_t>::const_iterator iter = _manual_modes.begin(); iter != _manual_modes.end(); ++iter)
        gradient(iter->first) = 0.0;
    }
    
    /** Maps random coefficients in the PCA space of the sample data to the corresponding \em shape and stores the squared Mahalanobis distance of the shape sample to the mean shape in \em squared_distance. The coefficients are sampled from a normal distribution with mean zero and covariance matrix <em>sample_radius * identity</em>. The manually set modes are not reflected in  \em squared_distance. */
    void random_shape_sample(float_t sample_radius, shape_t & out_shape, float_t & squared_distance) const
    {
      ublas::vector<float_t> coefficients(_covariance.size1());
      
      for(std::size_t i = 0; i < coefficients.size(); ++i)
      {
        if(_manual_modes.count(i) > 0)
          coefficients(i) = symmetric_uniform_distribution();
        else
          coefficients(i) = sample_radius * normal_distribution();
      }
          
      shape_sample(coefficients, out_shape, squared_distance);
    }
    
    /** Returns the mean shape. The user must call compute_statistics() before the mean shape can be queried. */
    const shape_t & mean_shape() const { return _mean_shape; }
    
    /** Sets \em out_shape to the shape on the <em>mode</em>-th principal component of the sample data at distance \em offset from the mean. Negative values for \em offset result in the opposite shape. */
    shape_t & mode_shape(std::size_t mode, float_t offset, shape_t & out_shape) const
    {
      float_t distance_from_mean;
      
      ublas::vector<float_t> coefficients(_covariance.size1());
      
      coefficients = ublas::scalar_vector<float_t>(coefficients.size(), 0.0);
        
      coefficients(mode) = offset;
      
      shape_sample(coefficients, out_shape, distance_from_mean);
      
      return out_shape;
    }
    
    /** Returns the dimension of the mean shape. */
    size_t dimension() const
    {
      return _mean_shape.dimension();
    }
    
  };
  
}


#endif
