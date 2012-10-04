// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef FEM_SIMPLEASSEMBLER_H
#define FEM_SIMPLEASSEMBLER_H

#include <fem/Grid.hpp>
#include <fem/FemKernel.hpp>
#include <fem/Assembler.hpp>
#include <fem/equation/SimpleEquationAdaptor.hpp>




namespace imaging
{
  /** \ingroup fem
      \brief Assembles the stiffness matrix and force vector of a FE problem implementing SimpleEquationInterface.
      
      The SimpleAssembler class provides functions to assemble the stiffness matrix and force vector for a given equation and a given grid.
      In contrast to Assembler, SimpleAssembler requires the equation to implement SimpleEquationInterface.
  */
  class SimpleAssembler
  {
    Assembler _assembler;
    
  public:

    /** Assembles the stiffness matrix and the force vector for \em simple_equation on \em grid. This is done in one big loop and thus faster than calling assemble_stiffness_matrix() and assemble_force_vector() separately. For performance reasons the type of \em simple_equation is a template parameter.
    
        The class \em simple_equation_t must implement all the functions defined in SimpleEquationInterface. It is also advised to derive \em simple_equation_t from SimpleEquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */
    template<class fem_types, class simple_equation_t>
    void assemble(const simple_equation_t & simple_equation, const Grid<fem_types> & grid,
                  ublas::compressed_matrix<float_t> & stiffness_matrix,
                  ublas::vector<float_t> & force_vector) const;
    
    /** Assembles the stiffness matrix for \c simple_equation on \c grid. If you want to assembly both, stiffness matrix and force vector, use assemble() to save computation time. For performance reasons the type of \em simple_equation is a template parameter.
    
        The class \em simple_equation_t must implement all the functions defined in EquationInterface. It is also advised to derive \em simple_equation_t from EquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */            
    template<class fem_types, class simple_equation_t>
    void assemble_stiffness_matrix(const simple_equation_t & simple_equation, const Grid<fem_types> & grid,
                  ublas::compressed_matrix<float_t> & stiffness_matrix) const;
                  
    /** Assembles the force vector for \c simple_equation on \c grid. If you want to assembly both, stiffness matrix and force vector, use assemble() to save computation time. For performance reasons the type of \em simple_equation is a template parameter.
    
        The class \em simple_equation_t must implement all the functions defined in EquationInterface. It is also advised to derive \em simple_equation_t from EquationInterface.
    
        The sparse matrix \em stiffness_matrix must be square and its size equal to the total number of nodes of the grid, as obtained from Grid::n_nodes(). It can greatly improve the performance of this function, if the matrix is pre-filled with zeros at positions where non-zero entries are expected. The position of the non-zeros entries depend on the geometry of the FE problem. Some grid construction functions provide an appropriate pre-filling of the stiffness matrix. 
    
        \sa Image2Grid, uniform_grid()
    */
    template<class fem_types, class simple_equation_t>
    void assemble_force_vector(const simple_equation_t & simple_equation, const Grid<fem_types> & grid,
                  ublas::vector<float_t> & force_vector) const;

  }
  ;

  template<class fem_types, class simple_equation_t>
  void SimpleAssembler::assemble(const simple_equation_t & simple_equation, const Grid<fem_types> & grid,
                                       ublas::compressed_matrix<float_t> & stiffness_matrix,
                                       ublas::vector<float_t> & force_vector) const
  {
    SimpleEquationAdaptor<simple_equation_t> adaptor(simple_equation);
    _assembler.assemble(adaptor, grid, stiffness_matrix, force_vector);
  }

  template<class fem_types, class simple_equation_t>
  void SimpleAssembler::assemble_stiffness_matrix(const simple_equation_t & simple_equation,
                                      const Grid<fem_types> & grid,
                                      ublas::compressed_matrix<float_t> & stiffness_matrix) const
  {
    SimpleEquationAdaptor<simple_equation_t> adaptor(simple_equation);
    _assembler.assemble_stiffness_matrix(adaptor, grid, stiffness_matrix);
  }

  template<class fem_types, class simple_equation_t>
  void SimpleAssembler::assemble_force_vector(const simple_equation_t & simple_equation,
                                      const Grid<fem_types> & grid,
                                      ublas::vector<float_t> & force_vector) const
  {
    SimpleEquationAdaptor<simple_equation_t> adaptor(simple_equation);
    _assembler.assemble_force_vector(adaptor, grid, force_vector);
  }
}


#endif
