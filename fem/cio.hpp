// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_CIO_H
#define FEM_CIO_H

#include <ostream>

#include "core/cio.hpp"
#include "fem/Element.hpp"
#include "fem/BoundaryElement.hpp"
#include "fem/Grid.hpp"

namespace imaging
{
  template<size_t N_VERTICES, size_t N_NODES>
  std::ostream & operator<<(std::ostream & out, const Element<N_VERTICES, N_NODES> & element)
  {
    for(int i = 0; i < N_VERTICES; ++i)
      out << element.get_vertex(i) << " ";
      
    out << "|";
    
    for(int i = 0; i < N_NODES; ++i)
      out << " " << element.get_node(i);
      
    out << " |";
    
    for(int i = 0; i < N_NODES; ++i)
      out << " " << element.get_node_status(i);
      
    out << std::endl;
    
    return out;
  }


  std::ostream & operator<<(std::ostream & out, const BoundaryElement & element)
  {
    out << element.get_parent_element() << " | " << element.get_parent_element_face() << std::endl;
    
    return out;
  }
  
 
  template<class fem_types>
  std::ostream & operator<<(std::ostream & out, const Grid<fem_types> & grid)
  {
    for(int i = 0; i < grid.get_n_vertices(); ++i)
      out << grid.get_vertex(i) << std::endl;
      
    for(int i = 0; i < grid.n_elements(); ++i)
      out << grid.get_element(i);
      
    for(int i = 0; i < grid.n_boundary_elements(); ++i)
      out << grid.get_boundary_element(i);
    
    return out;
  }

}


#endif
