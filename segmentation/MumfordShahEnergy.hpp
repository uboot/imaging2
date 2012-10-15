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

#ifndef ENERGY_MUMFORDSHAH_H
#define ENERGY_MUMFORDSHAH_H

#include <image/Image.hpp>
#include <image/ScalarImage.hpp>
#include <image/utilities.hpp>
#include <shape/BoundaryDiscretizer.hpp>
#include <shape/ShapeStatistics.hpp>
#include <shape/ShapeEnergyInterface.hpp>
#include <minimize/DifferentiableEnergyInterface.hpp>


namespace imaging 
{
  /** \ingroup segmentation
      \brief Simplified Mumford-Shah energy for arbitrary shape classes.
      
      This class template implements the simplified version of the Mumford-Shah functional for shapes of the type \em shape_t.
      The template parameter \em shape_t must implement ShapeInterface. 
      
      The simplified Mumford-Shah functional reads as follows:
      \f[
      I_\beta^{\textrm{SMS}}(\gamma) = \int_{\mathcal{I}(\gamma)}\big(u_1(\gamma) - f\big)^2 d
      x + \int_{\mathcal{O}(\gamma)} 
      \big(u_2(\gamma) - f\big)^2 d x + \beta |\gamma|_2^2\,,  
      \f]
      where
      \f[
      u_1(C) = \frac1{|\mathcal{I}(\gamma)|}\int_{\mathcal{I}(\gamma)}f d x
      \quad{\textrm{and}} \quad
      u_2(C) = \frac1{|\mathcal{O}(\gamma)|}\int_{\mathcal{O}(\gamma)}f d x\,.
      \f]
      
      Upon construction a MumfordShah object is initialized to an Image \f$u\f$ and an initial shape \f$\mu\f$.
      
      Then, the MumfordShah objects maps a vector \f$v\f$ to the corresponding energy as follows:
      \f[
      v \mapsto I_\beta^{\textrm{SMS}}\big(\textrm{Exp}_\mu(v)\big)\,. 
      \f]
  */
  template <class shape_t>
  class MumfordShahEnergy : public DifferentiableEnergyInterface, public ShapeEnergyInterface<shape_t>
  {
    const static std::size_t N = shape_t::SHAPE_DIMENSION;
    
    Image<N, ublas::fixed_vector<float_t, N> > _contrast_vector_field;
    Image<N, ublas::fixed_vector<float_t, N> > _squared_contrast_vector_field;
    Image<N, ublas::fixed_vector<float_t, N> > _volume_vector_field;
    Image<N, float_t> _image;
    float_t _beta;
    std::size_t _n_integration_points;
    
    ublas::vector<float_t> _current_argument;
    ublas::vector<float_t> _current_gradient;
    float_t _current_energy;
    shape_t _current_shape;
    shape_t _initial_shape;
    
    float_t _image_volume;
    float_t _image_contrast;
    float_t _squared_image_contrast;
    
    void compute_energy(float_t & u_1, float_t & u_2, float_t & energy);
       
  public:
  
    /** Constructs a MumfordShahEnergy energy object from an image accessor and shape statistics. */
    template <class const_accessor_t>
    MumfordShahEnergy(const const_accessor_t & image,
                const shape_t & initial_shape,
                float_t beta,
                std::size_t n_integration_points);
    
    ublas::vector<float_t> & current_argument() { return _current_argument; }
    float_t current_energy() const { return _current_energy; }
    std::size_t dimension() const { return _initial_shape.dimension(); }
    const shape_t & current_shape() const { return _current_shape; }

    const ublas::vector<float_t> & current_gradient() const { return _current_gradient; }
    
    void set_argument();
    void set_argument_with_gradient();
  };
  
  template <class shape_t>
  template <class const_accessor_t>
  MumfordShahEnergy<shape_t>::MumfordShahEnergy(const const_accessor_t & image,
                const shape_t & initial_shape,
                float_t beta,
                std::size_t n_integration_points) :
    _beta(beta),
    _initial_shape(initial_shape),
    _current_shape(initial_shape),
    _n_integration_points(n_integration_points),
    _current_argument(ublas::scalar_vector<float_t>(initial_shape.dimension(), 0.0)),
    _current_gradient(initial_shape.dimension()),
    _contrast_vector_field(image.size()),
    _squared_contrast_vector_field(image.size()),
    _volume_vector_field(image.size()),
    _image(image)
  {
    Image<N, float_t> squared_image(_image.size());
    ublas::fixed_vector<size_t, N> index;
    index.assign(0);
    do
    {
      squared_image[index] = square(_image[index]);
    }
    while(increment_index(_image.size(), index));
  
    compute_divergence_field(image, _contrast_vector_field);
    compute_divergence_field(squared_image, _squared_contrast_vector_field);
    compute_divergence_field(ScalarImage<N, float_t>(image.size(), 1.0), _volume_vector_field);
    
    _image_volume = 1.0;
    
    for(std::size_t i = 0; i < N; ++i)
      _image_volume *= float_t(image.size()(i));
      
    _image_contrast = 0.0;
    _squared_image_contrast = 0.0;
    
    index.assign(0);
    do
    {
      _image_contrast += _image[index];
      _squared_image_contrast += squared_image[index];
    }
    while(increment_index(_image.size(), index));
                                                        
    set_argument_with_gradient();
  }
  
  template <class shape_t>
  void MumfordShahEnergy<shape_t>::set_argument()
  {
    float_t u_1, u_2;
    
    compute_energy(u_1, u_2, _current_energy);
  }
  
  
  template <class shape_t>
  void MumfordShahEnergy<shape_t>::set_argument_with_gradient()
  {
    shape_t perturbed_shape;
    ublas::vector<float_t> perturbed_argument(dimension());
    ublas::vector<float_t> statistics_gradient(dimension());
    const float_t FINITE_H = 0.001;
    ublas::fixed_vector<float_t, shape_t::SHAPE_DIMENSION> point, normal;
    float_t u_1, u_2;
    
    perturbed_argument = _current_argument;
    _current_gradient = ublas::scalar_vector<float_t>(dimension(), 0.0);
    
    compute_energy(u_1, u_2, _current_energy);
    
    std::auto_ptr< BoundaryDiscretizer<shape_t::SHAPE_DIMENSION> > discretizer = _current_shape.boundary_discretizer(_n_integration_points);

    for(std::size_t k = 0; k < dimension(); ++k)
    {
      perturbed_argument(k) += FINITE_H;
      
      _initial_shape.exponential(perturbed_argument, perturbed_shape);
      
      std::auto_ptr< BoundaryDiscretizer<shape_t::SHAPE_DIMENSION> > perturbed_discretizer = perturbed_shape.boundary_discretizer(_n_integration_points);

      ublas::fixed_vector<float_t, shape_t::SHAPE_DIMENSION> perturbed_point, perturbed_normal, shape_derivative, shape_normal_derivative;
      float_t image_value;
  
      for(std::size_t j = 0; j < _n_integration_points; ++j)
      {
        perturbed_point = (*perturbed_discretizer)(j, perturbed_normal);
        point = (*discretizer)(j, normal);
        
        shape_derivative = 1.0 / FINITE_H * (perturbed_point - point);
        shape_normal_derivative = 1.0 / FINITE_H * (perturbed_normal - normal);
        
        // explicit cast of floating point values in "point" to pixel position!
        ublas::fixed_vector<size_t, N> pixel_position = ublas::fixed_vector<size_t, N>(point);   
        bool is_valid_pixel = true;
        
        for(std::size_t i = 0; i < N; ++i)
        {
          if( ! ( pixel_position(i) < _image.size()(i) ) )
          {
            is_valid_pixel = false;
            break;
          }
        }
        
        if(is_valid_pixel)
        {
        
          image_value = _image[pixel_position];
          
          _current_gradient(k) += ( square(u_1 - image_value) - square(u_2 - image_value) ) *
                                  inner_prod(normal, shape_derivative) * sqrt(norm_2(normal));
        }
        _current_gradient(k) += 2.0 * _beta * inner_prod(normal, shape_normal_derivative);
      }  
        
      perturbed_argument(k) -= FINITE_H;   
    }
  }
  
  template <class shape_t>
  void MumfordShahEnergy<shape_t>::compute_energy(float_t & u_1, float_t & u_2, float_t & energy)
  {
    _current_gradient = ublas::scalar_vector<float_t>(dimension(), 0.0);
    
    _initial_shape.exponential(_current_argument, _current_shape);
    
    float_t inner_contrast, outer_contrast;
    float_t inner_squared_contrast, outer_squared_contrast;
    float_t inner_volume, outer_volume;
    
    std::auto_ptr< BoundaryDiscretizer<shape_t::SHAPE_DIMENSION> > discretizer = _current_shape.boundary_discretizer(_n_integration_points);
    
    inner_contrast = discretizer->integrate_vector_field(_contrast_vector_field);
    inner_squared_contrast = discretizer->integrate_vector_field(_squared_contrast_vector_field);
    inner_volume = discretizer->integrate_vector_field(_volume_vector_field);
    
    outer_contrast = _image_contrast - inner_contrast;
    outer_squared_contrast = _squared_image_contrast - outer_squared_contrast;
    outer_volume = _image_volume - inner_volume;
    
    u_1 = inner_contrast / inner_volume;
    u_2 = outer_contrast / outer_volume;
    
    float_t shape_boundary_energy = 0.0;
    ublas::fixed_vector<float_t, shape_t::SHAPE_DIMENSION> point, normal;
    for(std::size_t j = 0; j < _n_integration_points; ++j)
    {
      point = (*discretizer)(j, normal);
      shape_boundary_energy += inner_prod(normal, normal);
    }
        
    energy = square(u_1) * inner_volume - 2 * u_1 * inner_contrast + inner_squared_contrast +
                      square(u_2) * outer_volume - 2 * u_2 * outer_contrast + outer_squared_contrast +
                      _beta * shape_boundary_energy;
  
  }

}


#endif
