// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef STATISTICAL_REGULARIZATION_H
#define STATISTICAL_REGULARIZATION_H

#include <shape/ShapeEnergyInterface.hpp>
#include <shape/ShapeStatistics.hpp>


namespace imaging
{
  /** \ingroup segmentation
      \brief Statistical regularization functional based on an arbitrary shape segmentation functional.
      
      This class implements the following functional:
      
      The Snakes functional reads as follows:
      \f[
      I_\alpha(p) = F(p) + \alpha d_{\mu, \Sigma}(\mu, p)^2 \,.
      \f]
      
      Upon construction a StatisticalRegularization object is initialized by
        - an arbitrary shape segmentation functional for shapes of type \em shape_t, \f$F\f$, 
        - the statistics which define the Mahalanobis distance \f$d_{\mu, \Sigma}\f$, and
        - the regularization parameter \f$\alpha\f$.
  */
  template <class shape_t>
  class StatisticalRegularization : public ShapeEnergyInterface<shape_t>
  {
    ShapeEnergyInterface<shape_t> & _energy;
    ShapeStatistics<shape_t> & _statistics;
    float_t _statistics_coefficient;
    
    ublas::vector<float_t> _current_argument;
    float_t _current_energy;
    
  public:
    /** Constructs a SnakesEnergy object. The arguments are explained in the class description. */
    StatisticalRegularization(ShapeEnergyInterface<shape_t> & energy, ShapeStatistics<shape_t> & statistics, float_t alpha)
      : _energy(energy), _statistics(statistics), _statistics_coefficient(alpha), _current_argument(ublas::scalar_vector<float_t>(energy.dimension(), 0.0))
    {
      if(_statistics.dimension() != _energy.dimension())
        throw Exception("Exception: dimensions of energy and statistics do not agree in StatisticalRegularization::StatisticalRegularization().");
      set_argument();
    }
    
    ublas::vector<float_t> & current_argument()
    {
      return _current_argument;
    }
    
    void set_argument()
    {
      float_t squared_distance;
      
      _statistics.shape_vector(_current_argument, _energy.current_argument(), squared_distance);
      _energy.set_argument();
      _current_energy = _energy.current_energy() + _statistics_coefficient * squared_distance;
      
//       std::cout << current_energy() << " " << _energy.current_energy() << " " << squared_distance << std::endl;
//       std::cout << _current_argument << std::endl;
//       std::cout << _energy.current_argument() << std::endl;
    }
    
    
    float_t current_energy() const
    {
      return _current_energy;
    }
    
    std::size_t dimension() const
    {
      return _energy.dimension();
    }
    
    const shape_t & current_shape() const
    {
      return _energy.current_shape();
    }
    
    

  };
}

#endif
