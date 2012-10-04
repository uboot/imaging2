// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef MREPMODEL_H
#define MREPMODEL_H

#include <shape/ShapeInterface.hpp>


namespace imaging
{
  /** \brief %MrepModel */
  template<class position_t, class atom_t, class connection_t>
  class MrepModel : public ShapeInterface
  {
  protected:
    position_t _position;

    std::vector<atom_t> _atoms;
    std::vector<connection_t> _connections;

  public:
    MrepModel() : _atoms(0), _connections(0) {}
    
    MrepModel(const position_t & position, size_t n_atoms = 0, size_t n_connections = 0) :
        _position(position), _atoms(n_atoms), _connections(n_connections)
    {}

    void resize(size_t n_atoms, size_t n_connections)
    {
      _atoms.resize(n_atoms);
      _connections.resize(n_connections);
    }

    const atom_t & atom(size_t i) const { return _atoms[i]; }
    void set_atom(size_t i, const atom_t & atom) { _atoms[i] = atom; }
    size_t n_atoms() const { return _atoms.size(); }

    const connection_t & connection(size_t i) const {return _connections[i]; }
    void set_connection(size_t i, const connection_t & connection) { _connections[i] = connection; }
    size_t n_connections() const { return _connections.size(); }

    const position_t & position() const {return _position; }
    void set_position(const position_t & position) { _position = position; }

    size_t start_atom(size_t atom_index) const
    {
      return _atoms[atom_index].start_atom();
    }

    size_t atom_connection(size_t atom_index) const
    {
      return _atoms[atom_index].connection();
    }


    void exponential(const ublas::vector<float_t> & vector, ShapeInterface & shape) const
    {
      MrepModel<position_t, atom_t, connection_t> * shape_ptr = dynamic_cast< MrepModel<position_t, atom_t, connection_t> *>(& shape);
      
      if(! shape_ptr)
        throw Exception("Exception: Wrong type of argument 'shape' in MrepModel::exponential().");

      shape_ptr->resize(n_atoms(), n_connections());

      ublas::vector<float_t>::const_iterator iter = vector.begin();
      _position.exponential(iter, shape_ptr->_position);

      for (size_t j = 0; j < n_atoms(); j++)
        _atoms[j].exponential(iter, shape_ptr->_atoms[j]);

      if (n_connections() > 0)
        _connections[0].exponential_without_rotation(iter, shape_ptr->_connections[0]);

      for (size_t j = 1; j < n_connections(); j++)
        _connections[j].exponential(iter, shape_ptr->_connections[j]);
    }


    void logarithm(const ShapeInterface & shape, ublas::vector<float_t> & vector) const
    {
      const MrepModel<position_t, atom_t, connection_t> * shape_ptr = dynamic_cast< const MrepModel<position_t, atom_t, connection_t> *>(& shape);
      
      if(! shape_ptr)
        throw Exception("Exception: Wrong type of argument 'shape' in MrepModel::logarithm().");
        
      if(n_atoms() != shape_ptr->n_atoms() || n_connections() != shape_ptr->n_connections())
        throw Exception("Exception: Wrong type of M-Rep in MrepModel::logarithm().");
      
      ublas::vector<float_t>::iterator iter = vector.begin();
      
      _position.logarithm(shape_ptr->_position, iter);

      for (size_t j = 0; j < n_atoms(); j++)
        _atoms[j].logarithm(shape_ptr->_atoms[j], iter);

      if (n_connections() > 0)
        _connections[0].logarithm_without_rotation(shape_ptr->_connections[0], iter);

      for (size_t j = 1; j < n_connections(); j++)
        _connections[j].logarithm(shape_ptr->_connections[j], iter);
    }

    size_t dimension() const
    {
      size_t dimension = 0;

      dimension += _position.dimension();

      for(size_t j = 0; j < n_atoms(); j++)
        dimension += _atoms[j].dimension();

      if(n_connections() > 0)
        dimension += _connections[0].dimension_without_rotation();

      for(size_t j = 1; j < n_connections(); j++)
        dimension += _connections[j].dimension();

      return dimension;
    }
  };


  /*template <class T, class P, class A, class C>
  std::ostream& operator<<( std::ostream & os, const MrepModel<T, P, A, C> & model)
  {
    os << "MrepModel";

    os << "\n";

    os << "  position: " << model.get_position() << "\n";
    os << "  atoms: " << model._atoms;
    os << "  connections: " << model._connections;

    os << "  atom table: " << model._atom_table;
    //     os << "  knot table: " << model._knot_table;

    os << std::endl;

    return os;
  }*/

}

#endif

