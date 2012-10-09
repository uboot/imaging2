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
