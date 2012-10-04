#include <shape/mrep/MrepModel2d.hpp>

#include <core/utilities.hpp>
#include <core/matrix_utilities.hpp>
#include <core/vector_utilities.hpp>

namespace imaging
{
  bool MrepModel2d::intersect( const ublas::fixed_vector<float_t, 2> & a_1, const ublas::fixed_vector<float_t, 2> & a_2,
                  const ublas::fixed_vector<float_t, 2> & b_1, const ublas::fixed_vector<float_t, 2> & b_2,
                  ublas::fixed_vector<float_t, 2> & intersection_point )
  {
    ublas::fixed_matrix<float_t, 2, 2> matrix;

    matrix(0,0) = a_2(0) - a_1(0);
    matrix(0,1) = - (b_2(0) - b_1(0));
    matrix(1,0) = a_2(1) - a_1(1);
    matrix(1,1) = - (b_2(1) - b_1(1));

    ublas::fixed_vector<float_t, 2> rhs(b_1(0) - a_1(0), b_1(1) - a_1(1));

    try
    {
      ublas::fixed_vector<float_t, 2> result;
      result = prod(inverse(matrix), rhs);

      if(result(0) >= 0.0 && result(0) < 1.0 && result(1) >= 0.0 && result(1) < 1.0 )
      {
        intersection_point = a_1 + result(0) * (a_2 - a_1);
        return true;
      }
      else
        return false;
    }
    catch(MathException)
    {
      return false;
    }
  }

  float_t MrepModel2d::tangent_factor(float_t angle)
  {
    return 2.5 * (1. + 1./PI * angle);
  }


  void MrepModel2d::init_boundary_coefficients(
    const std::vector< ublas::fixed_vector<float_t, 2> > & position_equations,
    const std::vector< ublas::fixed_vector<float_t, 2> > & tangent_equations,
    curve_t & spline_curve)
  {
    size_t interpolation_points = position_equations.size();

    /* Initialize the knot vector. Note that these are "knots" in the sense
       of splines, i.e. real values in the parameter domain of the spline 
       curve. In the MRep classes we often speak of the boundary postd::size_ts
       as knots. These are called "positions" here. */

    spline_curve.resize(SPLINE_ORDER, 2 * interpolation_points);

    for(size_t j = 0; j < interpolation_points; j++)
    {
      /* Every knot is inserted twice. I.e. the resulting spline curve is C^1 only
         On the other hand this results in more degrees of freedom for the 
         std::size_terpolation. */
      spline_curve.set_knot(j * 2, j);
      spline_curve.set_knot(j * 2 + 1, j);
    }

    spline_curve.set_knot(interpolation_points * 2, interpolation_points);

    /* Compute the spline coefficients from the position and tangent information.
       Note that we essentially solve the spline std::size_terpolation problem here. Because
       of the special structure of the spline we can explicitly compute the coefficients.
       Under more general circumstances we would have to solve a system of linear equations.
       This would happen, if we used the function BSplineCurve::std::size_terpolate(). */

    for(size_t j = 0; j < interpolation_points; j++)
    {
      spline_curve.set_coefficient(2*j, position_equations[j] -
                                   1.0/3.0 * tangent_equations[j]);

      spline_curve.set_coefficient(2*j + 1, position_equations[j] +
                                   1.0/3.0 * tangent_equations[j]);
    }
  }
  

  void MrepModel2d::compute_boundary(curve_t & spline_curve) const
  {
    // All of the following 3 vectors have exactly one entry per atom.

    /* A list of the absolute positions of all atoms (in the same order as in
       MRepModel::_atoms). */
    std::vector< ublas::fixed_vector<float_t, 2> > atom_positions(n_atoms());

    /* Here for every atom a list of its neighbor atoms is stored.
       I.e. in atom_targets[i] you can find a vector of the indices of atoms,
       which are connected to the i-th atom (by means of parts of the medial
       skeleton). */
    std::vector< std::vector<size_t> > atom_targets(n_atoms());

    /* If one atom is connected to another, we can define the direction of the
       connecting vector. These direction are stored in atom_target_angles.
       I.e. atom_targets_angles[i][j] is the direction of the connection connecting
       the i-th atom to the atom_targets[i][j]-th atom. */
    std::vector< std::vector<float_t> > atom_target_angles(n_atoms());

    // Write exactly this information std::size_to the vectors:
    for(size_t i = 1; i < n_atoms(); i++)
    {
      size_t next_atom;

      next_atom = start_atom(i);

      if(next_atom != i)
      {
        ublas::fixed_vector<float_t, 2> this_position = atom_center(i);
        ublas::fixed_vector<float_t, 2> next_position = atom_center(next_atom);

        atom_positions[i] = atom_center(i);

        atom_targets[i].push_back(next_atom);
        atom_targets[next_atom].push_back(i);

        atom_target_angles[i].push_back(angle(next_position - this_position));
        atom_target_angles[next_atom].push_back(angle(this_position - next_position));

      }
    }

    /* All of the following 2 vectors will have exactly one entry for every contact of an atom
       with the boundary. According to the idea of M-Reps there should be one
       contact for end (boundary) atoms and 2 for inner atoms. */

    /* The entry atom_list[i] corresponds to the atom index of the i-th contact
       between boundary and atoms. Basically, every contact will result in a spline postd::size_t. */
    std::vector<size_t> atom_list;

    /* Every contact corresponds to a tangent at the touching atom. Its direction
       is stored here. */
    std::vector<float_t> tangent_angle_list;

    /* In the following loop the vectors defined above are filled with the appropriate
       values. Note that here we connect atom boundary postd::size_ts by tangents only, i.e. 
       every connecting segment is line segment. In some cases these line segments
       will intersect before the touch their target atom. These cases will require
       some post-processing, when computing the final spline postd::size_ts. */

    // Start at the atom with index 0.
    size_t atom_index = 0;
    float_t incoming_angle(atom_target_angles[0][0]);
    size_t incoming_atom_index = atom_targets[0][0];
    size_t first_incoming_atom = incoming_atom_index;

    do
    {
      atom_list.push_back(atom_index);

      // Compute the index of the next atom.

      // If their are more than one atoms connected, we have to choose the right one.
      if(atom_target_angles[atom_index].size() > 1)
      {
        float_t result_angle;

        size_t target_index;

        if(atom_targets[atom_index][0] != incoming_atom_index)
          target_index = 0;
        else
          target_index = 1;

        size_t result_target_index = target_index;

        result_angle = clockwise_difference(incoming_angle, atom_target_angles[atom_index][target_index]);

        target_index++;

        // The right one is the one, s.t. the angle between incoming and outgoing tangent is smallest.
        for(; target_index < atom_targets[atom_index].size(); target_index++)
        {
          if(atom_targets[atom_index][target_index] != incoming_atom_index)
          {

            float_t angle = clockwise_difference(incoming_angle, atom_target_angles[atom_index][target_index]);

            if(angle < result_angle)
            {
              result_target_index = target_index;
              result_angle = angle;
            }
          }
        }

        incoming_angle = 2 * PI - atom_target_angles[atom_index][result_target_index];

        // Store the index of the right atom.
        incoming_atom_index = atom_index;
        atom_index = atom_targets[atom_index][result_target_index];
      }
      else
      {
        // Things are easy if there is only one atom connected.
        incoming_angle = atom_target_angles[atom_index][0] + PI;
        if(incoming_angle > 2 * PI) incoming_angle -= 2 * PI;
        incoming_atom_index = atom_index;
        atom_index = atom_targets[atom_index][0];
      }

      // Also store the direction of the tangent of this boundary contact.
      tangent_angle_list.push_back(compute_tangent_angle(incoming_atom_index, atom_index));
    }
    while(atom_index != 0 || incoming_atom_index != first_incoming_atom);


    // We are finally ready to compute the final values of the spline curve points and its tangents.
    std::vector< ublas::fixed_vector<float_t, 2> > knot_position_list;
    std::vector< ublas::fixed_vector<float_t, 2> > knot_tangent_list;

    // Again we have to start somewhere (1st contact of the spline curve).
    float_t angle_1 = tangent_angle_list[tangent_angle_list.size() - 1];
    float_t angle_2 = tangent_angle_list[0];

    ublas::fixed_vector<float_t, 2> in_point(atom_center(atom_list[atom_list.size() - 1]) +
                       polar2cartesian(atom_radius(atom_list[atom_list.size() - 1]), angle_1));

    ublas::fixed_vector<float_t, 2> current_point_1(atom_center(atom_list[0]) +
                              polar2cartesian(atom_radius(atom_list[0]), angle_1));

    ublas::fixed_vector<float_t, 2> current_point_2(atom_center(atom_list[0]) +
                              polar2cartesian(atom_radius(atom_list[0]), angle_2));

    ublas::fixed_vector<float_t, 2> out_point;

    // Start the loop.
    for(size_t j = 0; j < atom_list.size(); j++)
    {
      /* We are now at the j-th contact of the spline curve to with atom. The index
         of this atom is atom_list[j]. The position of the point, where the previous
         atom was touched, is stored in in_postd::size_t. We compute the position of the
         point, where the next atom is touched, and store it in out_point. 
         The current atom is touched at two points:
          (1) by the tangent from the previous atom (current_point_2), 
          (2) by the tangent to the next atom (current_point_1). */
      out_point = atom_center(atom_list[(j + 1) % atom_list.size()]) +
                  polar2cartesian(atom_radius(atom_list[(j + 1) % atom_list.size()]), angle_2);

      /* Thus, we get two line segments:
           (1) from the previous atom to the current one,
           (2) from the current one to the next one.
         Do the intersect? */

      ublas::fixed_vector<float_t, 2> intersection_point;

      if( ! intersect(current_point_1, in_point, current_point_2, out_point, intersection_point))
      {
        // They do not.
        float_t angle_difference = clockwise_difference(angle_1, angle_2);

        float_t radial_vector_angle = angle_1 - angle_difference / 2.0;

        knot_position_list.push_back(atom_center(atom_list[j]) +
                                     polar2cartesian(atom_radius(atom_list[j]), radial_vector_angle));

        knot_tangent_list.push_back(atom_radius(atom_list[j]) *
                                    polar2cartesian(tangent_factor(angle_difference), radial_vector_angle - PI / 2.0));
      }
      else
      {
        // They intersect in intersection_point.
        float_t angle_difference = counter_clockwise_difference(angle_1, angle_2);

        float_t radial_vector_angle = angle_1 + angle_difference / 2.0;

        knot_position_list.push_back(atom_center(atom_list[j]) +
                                     polar2cartesian(atom_radius(atom_list[j]), radial_vector_angle));

        knot_tangent_list.push_back(atom_radius(atom_list[j]) *
                                    polar2cartesian(tangent_factor(- angle_difference), radial_vector_angle - PI / 2.0));
      }


      // Prepare to touch the next atom!
      angle_1 = angle_2;
      angle_2 = tangent_angle_list[(j + 1) % tangent_angle_list.size()];

      in_point = current_point_2;
      current_point_1 = out_point;
      current_point_2 = atom_center(atom_list[(j + 1) % atom_list.size()]) +
                        polar2cartesian(atom_radius(atom_list[(j + 1) % atom_list.size()]), angle_2);
    }

    // Finally compute the spline coefficients from the boundary postd::size_t and tangent information.
    init_boundary_coefficients(knot_position_list, knot_tangent_list, spline_curve);
  }

  
  MrepModel2d::Discretizer::Discretizer(const MrepModel2d & model, size_t n_points) : img::BoundaryDiscretizer<2>(n_points)
  {
    model.compute_boundary(_boundary);
    _step_size = ( _boundary.last_knot() - _boundary.first_knot() ) / float_t(_n_points);
  }

  void MrepModel2d::Discretizer::evaluate(size_t i, ublas::fixed_vector<float_t, 2> & point, ublas::fixed_vector<float_t, 2> & normal, float_t & curvature) const
  {
    ublas::fixed_vector<float_t, 2> tangent, second_derivative;
    point = _boundary(i * _step_size, tangent, second_derivative);

    tangent *= _step_size;

    normal(0) = - tangent(1);
    normal(1) = tangent(0);
    
    curvature = boundary_discretizer_impl::compute_curve_curvature(tangent, second_derivative);
  }
  
  std::auto_ptr< BoundaryDiscretizer<2> > MrepModel2d::boundary_discretizer(size_t n_points) const
  {
    return std::auto_ptr< BoundaryDiscretizer<2> >(new Discretizer(*this, n_points));
  }

}

