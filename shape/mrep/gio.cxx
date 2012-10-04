#include <shape/mrep/gio.hpp>

#include <polytope/gio.hpp>
#include <spline/gio.hpp>

namespace imaging
{
  GraphicsInterface & operator<<(GraphicsInterface & out, const MrepModel2d & model)
  {
    MrepModel2d::curve_t boundary;

    model.compute_boundary(boundary);

    GraphicsInterface::StreamStatus stream_status_backup = out.get_stream_status();

    out << GraphicsInterface::set_line_width(2.0) << GraphicsInterface::set_color(Color::RED);
    out << boundary;

    out << GraphicsInterface::set_color(Color::YELLOW);
    for(size_t i = 0; i < model.n_atoms(); ++i)
      out.vertex(model.atom_center(i));

    out << GraphicsInterface::set_line_width(1.0);
    out << GraphicsInterface::set_color(Color::CYAN);
    out << GraphicsInterface::offset_z_layer(-1);

    for(size_t i = 0; i < model.n_atoms(); ++i)
      out.circle(model.atom_center(i), model.atom_radius(i));

    out << GraphicsInterface::set_color(Color::BLUE);
    out << GraphicsInterface::offset_z_layer(-1);

    std::vector< ublas::fixed_vector<float_t, 2> > line_segment(2);
    for(size_t i = 0; i < model.n_atoms(); ++i)
    {
      line_segment[0] = model.atom_center(i);
      line_segment[1] = model.atom_center(model.atom(i).start_atom());

      out.polyline(line_segment);
    }

    out << GraphicsInterface::offset_z_layer(+2);
    out << stream_status_backup;

    return out;
  }
  
  
  GraphicsInterface & operator<<(GraphicsInterface & out, const PolygonModel2d & model)
  {
    Polygon polygon;
    
    model.compute_boundary(5, polygon);

    GraphicsInterface::StreamStatus stream_status_backup = out.get_stream_status();
    
    out << GraphicsInterface::set_line_width(2.0) << GraphicsInterface::set_color(Color::RED);
    out << polygon;
    
    out << GraphicsInterface::set_color(Color::YELLOW);
    for(size_t i = 0; i < model.n_atoms(); ++i)
      out.vertex(model.atom_center(i));

    out << GraphicsInterface::set_line_width(1.0);
    out << GraphicsInterface::set_color(Color::CYAN);
    out << GraphicsInterface::offset_z_layer(-1);

    for(size_t i = 0; i < model.n_atoms(); ++i)
      out.circle(model.atom_center(i), model.atom_radius(i));

    out << GraphicsInterface::set_color(Color::BLUE);
    out << GraphicsInterface::offset_z_layer(-1);

    std::vector< ublas::fixed_vector<float_t, 2> > line_segment(2);
    for(size_t i = 0; i < model.n_atoms(); ++i)
    {
      line_segment[0] = model.atom_center(i);
      line_segment[1] = model.atom_center(model.atom(i).start_atom());

      out.polyline(line_segment);
    }
    
//     std::auto_ptr< BoundaryDiscretizer<2> > discretizer = model.boundary_discretizer(50);
//     for(size_t i = 0; i < discretizer->n_points(); ++i)
//     {
//       ublas::fixed_vector<float_t, 2> point, normal;
//       point = (*discretizer)(i, normal);
//       out.arrow(point, normal);
//     }

    out << GraphicsInterface::offset_z_layer(+2);
    out << stream_status_backup;
    
    return out;
  }
}

