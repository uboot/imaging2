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

#ifndef SNAKESENERGY_H
#define SNAKESENERGY_H

#include <shape/ShapeEnergyInterface.hpp>
#include <image/Image.hpp>
#include <shape/ShapeStatistics.hpp>
#include <shape/BoundaryDiscretizer.hpp>

#include <set>


namespace imaging
{
  /** \ingroup segmentation
      \brief Snakes energy for arbitrary shape classes.
      
      This class template implements the Snakes functional for shapes of the type \em shape_t.
      The template parameter \em shape_t must implement ShapeInterface. 
      
      The Snakes functional reads as follows:
      \f[
      I_\beta^{\textrm{Snakes}}(\gamma) = - \int_\gamma u(\gamma) d \tau + 
      \beta \textrm{Length}(\gamma)\,.
      \f]
      
      Usually \f$u\f$ is chosen as 
      \f[
      u = | \nabla f |\,,
      \f]
      for an image \f$f\f$. Thus, \f$u\f$ contains the edges of \f$f\f$.
      
      Upon construction a SnakesEnergy object is initialized to an Image \f$u\f$ and an initial shape \f$\mu\f$.
      
      Then, the SnakesEnergy objects maps a vector \f$v\f$ to the corresponding energy as follows:
      \f[
      v \mapsto I_\beta^{\textrm{Snakes}}\big(\textrm{Exp}_\mu(v)\big)\,. 
      \f]
  */
  template <class shape_t>
  class SnakesEnergy : public ShapeEnergyInterface<shape_t>
  {
    const static std::size_t SHAPE_DIMENSION = shape_t::SHAPE_DIMENSION;
    
    Image<SHAPE_DIMENSION, float_t> _edge_map;
    float_t _beta;
    std::size_t _n_integration_points;
    
    shape_t _initial_shape;
    ublas::vector<float_t> _current_argument;
    shape_t _current_shape;
    float_t _current_energy;
    
  public:
    /** Constructs a SnakesEnergy object. The first three arguments are explained in the class description. The parameter \em n_integration_points refers to the number of integration points when integrating along the shape boundary. */
    template <class const_accessort_t>
    SnakesEnergy(const const_accessort_t & edge_map,
                 const shape_t & initial_shape,
                 float_t beta,
                 std::size_t n_integration_points) :
    _edge_map(edge_map),
    _initial_shape(initial_shape), 
    _beta(beta),
    _n_integration_points(n_integration_points),
    _current_argument(ublas::scalar_vector<float_t>(initial_shape.dimension(), 0.0)),
    _current_shape(initial_shape),
    _current_energy(0.0)
    { }
    
    ublas::vector<float_t> & current_argument()
    {
      return _current_argument;
    }
    
    const shape_t & current_shape() const
    {
      return _current_shape;
    }
    
    float_t current_energy() const
    {
      return _current_energy;
    }
    
    std::size_t dimension() const
    {
      return _initial_shape.dimension();
    }
    
    void set_argument()
    {
      float_t boundary_area, edge_energy;
      std::set<size_t> overlap_indices;
      
      _initial_shape.exponential(_current_argument, _current_shape);
      
      std::auto_ptr< BoundaryDiscretizer<shape_t::SHAPE_DIMENSION> > edge_map_discretizer = _current_shape.boundary_discretizer(_n_integration_points);

      boundary_area = edge_map_discretizer->compute_boundary_area();
      edge_energy = edge_map_discretizer->integrate(_edge_map);
      
      _current_energy = - edge_energy + _beta * boundary_area;
                    
    }           
  };
}

#endif
