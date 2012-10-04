// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_MREP_MREPMODEL2D_H
#define SHAPE_MREP_MREPMODEL2D_H

#include <shape/mrep/MrepSkeleton2d.hpp>
#include <spline/PeriodicBspline.hpp>
#include <shape/DiscretizableShapeInterface.hpp>


namespace imaging
{
  /** \ingroup mrep
      \brief Medial axis parametrization of planar shapes.
  */
  class MrepModel2d : public MrepSkeleton2d, public DiscretizableShapeInterface<2>
  {
  public:
    typedef PeriodicBspline< ublas::fixed_vector<float_t, 2> > curve_t;

  private:
    class Discretizer;
    
    static const size_t SPLINE_ORDER = 4;
    
    static void init_boundary_coefficients(
      const std::vector< ublas::fixed_vector<float_t, 2> > & position_equations,
      const std::vector< ublas::fixed_vector<float_t, 2> > & tangent_equations,
      curve_t & spline_curve);
    static float_t tangent_factor(float_t angle);
    static bool intersect( const ublas::fixed_vector<float_t, 2> & a_1, const ublas::fixed_vector<float_t, 2> & a_2,
                    const ublas::fixed_vector<float_t, 2> & b_1, const ublas::fixed_vector<float_t, 2> & b_2,
                    ublas::fixed_vector<float_t, 2> & intersection_point );



  public:
    MrepModel2d() : MrepSkeleton2d() {}
    MrepModel2d(const Position2d & position,
                size_t n_atoms = 0, size_t n_connections = 0) :
      MrepSkeleton2d(position, n_atoms, n_connections) {}
      
    virtual std::auto_ptr< BoundaryDiscretizer<2> > boundary_discretizer(size_t n_points) const;

    //! Computes the boundary spline curve of the model.
    void compute_boundary(curve_t & spline_curve) const;
  };
    
  /** \cond */
  class MrepModel2d::Discretizer : public BoundaryDiscretizer<2>
  {
    PeriodicBspline< ublas::fixed_vector<float_t, 2> > _boundary;
    float_t _step_size;

  public:
    
    Discretizer(const MrepModel2d & model, size_t n_points);

    void evaluate(size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & point, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const;
  };
  /** \endcond */
  
}

#endif
