#include <shape/mrep/MrepSkeleton2d.hpp>

#include <core/utilities.hpp>
#include <core/matrix_utilities.hpp>
#include <core/vector_utilities.hpp>

namespace imaging
{
  float_t MrepSkeleton2d::compute_tangent_angle(size_t atom_1, size_t atom_2) const
  {
    ublas::fixed_vector<float_t, 2> center_1 = atom_center(atom_1);
    ublas::fixed_vector<float_t, 2> center_2 = atom_center(atom_2);
    float_t radius_1 = atom_radius(atom_1);
    float_t radius_2 = atom_radius(atom_2);

    float_t result;

    ublas::fixed_vector<float_t, 2> d;
    float_t alpha, cos_alpha;

    d = center_2 - center_1;

    cos_alpha = (radius_2 - radius_1)/radius(d);
    if(cos_alpha > 1.0) cos_alpha = 1.0;
    else if(cos_alpha < -1.0) cos_alpha = -1.0;

    alpha = PI - acos(cos_alpha);

    result = angle(d) + alpha;

    if(result < 0)
      result += 2 * PI;

    return result;
  }
  
  
  float_t MrepSkeleton2d::atom_rotation(size_t atom_index) const
  {
    float_t rotation(0.0);
    size_t atom = atom_index;

    while(atom != 0)
    {
      size_t connection_index = atom_connection(atom);
      rotation += connection(connection_index).rotation();

      atom = start_atom(atom);
    }

    return position().rotation() + rotation;
  }


  ublas::fixed_vector<float_t, 2> MrepSkeleton2d::atom_center(size_t atom_index) const
  {
    ublas::fixed_vector<float_t, 2> center(0, 0);
    size_t atom = atom_index;

    while(atom != 0)
    {
      float_t start_atom_rotation = atom_rotation(start_atom(atom));

      float_t min_length = fabs(atom_radius(start_atom(atom)) - atom_radius(atom));
      float_t length =  min_length + connection(atom_connection(atom)).radius();

      ublas::fixed_vector<float_t, 2> temp;
      ublas::fixed_matrix<float_t, 2, 2> rotation;
 
      if(atom_connection(atom) != 0)
        rotation = rotation_matrix(start_atom_rotation + connection(atom_connection(atom)).rotation());
      else
        rotation = rotation_matrix(start_atom_rotation);
        
      
      center += length * prod(rotation, ublas::fixed_vector<float_t, 2>(1.0, 0.0));

      atom = start_atom(atom);
    }

    return position().center() + center;
  } 
  
  void MrepSkeleton2d::set_atom_center(size_t atom_index, const ublas::fixed_vector<float_t, 2> & center)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > atom_centers(n_atoms());
    std::vector<float_t> atom_radii(n_atoms());
    
    get_geometry(atom_centers, atom_radii);
    
    atom_centers[atom_index] = center;
    
    set_geometry(atom_centers, atom_radii);
  }
  
  void MrepSkeleton2d::set_atom_radius(size_t atom_index, float_t radius)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > atom_centers;
    std::vector<float_t> atom_radii;
    
    get_geometry(atom_centers, atom_radii);
    
    atom_radii[atom_index] = radius;
    
    set_geometry(atom_centers, atom_radii);
  }

  void MrepSkeleton2d::get_geometry(std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers,
    std::vector<float_t> & atom_radii) const
  {
    atom_centers.resize(n_atoms());
    atom_radii.resize(n_atoms());
    
    for(size_t i = 0; i < n_atoms(); ++i)
    {
      atom_centers[i] = atom_center(i);
      atom_radii[i] = atom_radius(i);
    }
  }
  
  void MrepSkeleton2d::set_geometry(const std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers,
    const std::vector<float_t> & atom_radii)
  {
    for(size_t i = 0; i < n_atoms(); ++i)
      _atoms[i].set_radius(atom_radii[i]);
    
    size_t atom_index;
    float_t rotation;
    
    for(atom_index = 1; atom_index < n_atoms(); ++atom_index)
      if(_atoms[atom_index].connection() == 0)
      {
        ublas::fixed_vector<float_t, 2> connection_vector =
          atom_centers[atom_index] - atom_centers[_atoms[atom_index].start_atom()];
        
        rotation = angle(connection_vector);
 
        set_position(Position2d(atom_centers[0], rotation));
        
        // _connections[0].set_radius(max(norm_2(connection_vector) - min_length, 0.0));
        // _connections[0].set_rotation(0.0);
      }
      
    set_subtree_geometry(0, rotation, atom_centers);
  }
  
  void MrepSkeleton2d::set_subtree_geometry(size_t atom_index, float_t rotation, const std::vector< ublas::fixed_vector<float_t, 2> > & atom_centers)
  {
    for(size_t i = 1; i < n_atoms(); ++i)
    {
      if(_atoms[i].start_atom() == atom_index)
      {
        ublas::fixed_vector<float_t, 2> connection_vector =
          atom_centers[i] - atom_centers[atom_index];
        
        float_t new_rotation = angle(connection_vector);
        float_t min_length = fabs(atom_radius(i) - atom_radius(atom_index));
        
        _connections[_atoms[i].connection()].set_radius(max(norm_2(connection_vector) - min_length, 0.0));
        _connections[_atoms[i].connection()].set_rotation(counter_clockwise_difference(rotation, new_rotation));
        
        set_subtree_geometry(i, new_rotation, atom_centers);
      }
    }
  }
}

