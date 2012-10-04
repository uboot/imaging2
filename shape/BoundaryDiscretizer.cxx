#include <shape/BoundaryDiscretizer.hpp>


extern "C" {
double inter(const double * a, int na, const double * b, int nb);
}


namespace imaging
{
  float_t compute_intersection_volume(const BoundaryDiscretizer<2> & shape_a, const BoundaryDiscretizer<2> & shape_b)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > points_a(shape_a.n_points());
    std::vector< ublas::fixed_vector<float_t, 2> > points_b(shape_b.n_points());
    
    if(shape_a.n_points() == 0 || shape_b.n_points() == 0)
      return 0.0;
      
    points_a[0] = shape_a(0);
    points_b[0] = shape_b(0);
    
    Box<2> bbox_a(points_a[0], points_a[0]);
    Box<2> bbox_b(points_b[0], points_b[0]);
    
    for(size_t i = 1; i < shape_a.n_points(); ++i)
    {
      points_a[i] = shape_a(i);
      bbox_a.add_point(points_a[i]);
    }
    
    for(size_t i = 1; i < shape_b.n_points(); ++i)
    {
      points_b[i] = shape_b(i);
      bbox_b.add_point(points_b[i]);
    }
    
    if(! do_intersect(bbox_a, bbox_b))
      return 0.0; 

    float_t result = inter(&(points_a[0](0)), shape_a.n_points(), &(points_b[0](0)), shape_b.n_points());
    
    return result;
  }
  
  namespace boundary_discretizer_impl
  {
    float_t compute_curve_curvature(ublas::fixed_vector<float_t, 2> & first_derivative, ublas::fixed_vector<float_t, 2> & second_derivative)
    {
      return abs(first_derivative(0) * second_derivative(1) - first_derivative(1) * second_derivative(0)) / 
             pow(norm_2(first_derivative), 3.0);
    }
  }
}

