// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_GIO_H
#define FEM_GIO_H

#include <graphics/GraphicsInterface.hpp>
#include <fem/Grid.hpp>

namespace imaging
{
  template <class fem_types>
  GraphicsInterface & operator<<(GraphicsInterface & out, const Grid<fem_types> & grid)
  {
    GraphicsInterface::StreamStatus stream_status_backup = out.get_stream_status();

    out << GraphicsInterface::set_color(Color::RED);
    
    for(int i = 0; i < grid.n_elements(); ++i)
    {                     
      std::vector< ublas::fixed_vector<float_t, fem_types::data_dimension> > element_vertices(fem_types::transform_t::n_element_vertices);
      
      for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
        element_vertices[j] = grid.vertex(i, j);
      
      out.polygon(element_vertices);
    
      out << GraphicsInterface::offset_z_layer(2);
      
      for(size_t j = 0; j < fem_types::shape_function_t::n_element_nodes; ++j)
      {
        if(grid.is_boundary_node(i, j))
        {
          out << GraphicsInterface::set_color(Color::BLUE);
//           out.vertex(grid.vertex(i, j));
          out.arrow(grid.vertex(i, j), grid.boundary_normal(i, j));
          out << GraphicsInterface::set_color(Color::RED);
        }
      }
    
      out << GraphicsInterface::offset_z_layer(-2);
    }
    
    out << GraphicsInterface::offset_z_layer(1);
      
    out << GraphicsInterface::set_color(Color::BLUE);
    
    for(int i = 0; i < grid.n_boundary_elements(); ++i)
    {
      size_t parent_element = grid.parent_element(i);
      size_t element_face = grid.parent_element_face(i);
      ublas::fixed_vector<size_t, 2> local_vertex_indices(fem_types::transform_t::face_vertex(element_face, 0),
                                                  fem_types::transform_t::face_vertex(element_face, 1));
      out.line_segment(grid.vertex(parent_element, local_vertex_indices(0)),
                       grid.vertex(parent_element, local_vertex_indices(1)));
    }
    
    out << GraphicsInterface::offset_z_layer(-1);

    out << stream_status_backup;
    
    return out;
  }

}


#endif
