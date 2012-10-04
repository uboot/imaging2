// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef SHAPE_MREP_POLYGONMODEL2D_H
#define SHAPE_MREP_POLYGONMODEL2D_H

#include <shape/mrep/MrepSkeleton2d.hpp>
#include <spline/PeriodicBspline.hpp>
#include <shape/DiscretizableShapeInterface.hpp>
#include <polytope/Polygon.hpp>


namespace imaging
{
  /** \ingroup mrep
      \brief Medial axis parametrization of planar shapes.
  */
  class PolygonModel2d : public MrepSkeleton2d, public DiscretizableShapeInterface<2>
  {
  public:
    typedef PeriodicBspline< ublas::fixed_vector<float_t, 2> > curve_t;

  private:
    class Discretizer;
    
    void compute_limb_polygon(size_t atom_1, size_t atom_2, size_t n_sector_discretization_points, Polygon & polygon) const;

  public:
    PolygonModel2d() : MrepSkeleton2d() {}
    PolygonModel2d(const Position2d & position,
                size_t n_atoms = 0, size_t n_connections = 0) :
      MrepSkeleton2d(position, n_atoms, n_connections) {}
      
    virtual std::auto_ptr< BoundaryDiscretizer<2> > boundary_discretizer(size_t n_points) const;
      
    //! Computes the boundary polygon of the model.
    void compute_boundary(size_t n_sector_discretization_points, Polygon & polygon) const;
  };
  
  /** \cond */
  class PolygonModel2d::Discretizer : public BoundaryDiscretizer<2>
  {
    float_t _step_size;
    std::vector< ublas::fixed_vector<float_t, 2> > _points;
    std::vector< ublas::fixed_vector<float_t, 2> > _normals;
    
    static const size_t N_SECTOR_DISCRETIZATION_POINTS = 5;

  public:
    
    Discretizer(const PolygonModel2d & model, size_t n_points);

    void evaluate(size_t i, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & point, ublas::fixed_vector<float_t, SHAPE_DIMENSION> & normal, float_t & curvature) const;
  };  
  /** \endcond */

}

#endif
