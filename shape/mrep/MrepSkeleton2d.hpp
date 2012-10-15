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

#ifndef SHAPE_MREP_MREPSKELETON2D_H
#define SHAPE_MREP_MREPSKELETON2D_H

#include <shape/mrep/MrepModel.hpp>
#include <shape/mrep/Position2d.hpp>
#include <shape/mrep/MrepAtom.hpp>
#include <shape/mrep/MrepConnection.hpp>
#include <shape/mrep/MrepModel.hpp>


namespace imaging
{
  /** \ingroup mrep
      \brief Medial axis parametrization of planar shapes.
  */
  class MrepSkeleton2d : public MrepModel< Position2d, MrepAtom, MrepConnection<2> >
  {
    float_t atom_rotation(size_t atom_index) const;
      
    void set_subtree_geometry(size_t atom_index, float_t rotation, const std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers);

  public:
    MrepSkeleton2d() : MrepModel< Position2d, MrepAtom, MrepConnection<2> >() {}
    MrepSkeleton2d(const Position2d & position,
                size_t n_atoms = 0, size_t n_connections = 0) :
      MrepModel< Position2d, MrepAtom, MrepConnection<2> >(position, n_atoms, n_connections) {}
    
    //! Get the global coordinates of the center of the atom \c atom_index.
    ublas::fixed_vector<float_t, 2> atom_center(size_t atom_index) const;
    
    //! Set the global coordinates of the center of the atom \c atom_index.
    void set_atom_center(size_t atom_index, const ublas::fixed_vector<float_t, 2> & center);
    
    //! Set the radius the atom \c atom_index.
    void set_atom_radius(size_t atom_index, float_t radius);

    void get_geometry(std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers,
      std::vector<float_t> & atom_radii) const;
    
    void set_geometry(const std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers,
      const std::vector<float_t> & atom_radii);
    
    //! Get the radius of the atom \c atom_index.
    float_t atom_radius(size_t atom_index) const
    {
      return _atoms[atom_index].radius();
    }
    
    float_t compute_tangent_angle(size_t atom_1, size_t atom_2) const;
  };

}

#endif
