#include <shape/mrep/xmlio.hpp>

namespace imaging
{
  const std::string xml_handler< Position2d >::element_name = "position";
  const std::string xml_handler<MrepAtom>::element_name = "atom";
  const std::string xml_handler< MrepConnection<2> >::element_name = "connection";
  const std::string xml_handler<MrepSkeleton2d>::element_name = "mrep_skeleton_2d";
  const std::string xml_handler<MrepModel2d>::element_name = "mrep_model_2d";
  const std::string xml_handler<PolygonModel2d>::element_name = "polygon_model_2d";

  void xml_handler< Position2d >::read_object(XmlReader & in, Position2d & object) const
  {
    ublas::fixed_vector<float_t, 2> center;
    float_t rotation;

    in >> XmlReader::attribute("center") >> center;
    in >> XmlReader::attribute("rotation") >> rotation;

    object.assign(center, rotation);
  }

  void xml_handler< Position2d >::write_object(const Position2d & object, XmlWriter & out) const
  {
    out << XmlWriter::attribute("center") << object.center();
    out << XmlWriter::attribute("rotation") << object.rotation();
  }

  void xml_handler<MrepAtom>::read_object(XmlReader & in, MrepAtom & object) const
  {
    float_t radius;
    std::size_t start_atom;
    std::size_t connection;

    in >> XmlReader::attribute("r") >> radius;
    in >> XmlReader::attribute("start_atom") >> start_atom;
    in >> XmlReader::attribute("connection") >> connection;

    object.assign(radius, start_atom, connection);
  }

  void xml_handler<MrepAtom>::write_object(const MrepAtom & object, XmlWriter & out) const
  {
    out << XmlWriter::attribute("r") << object.radius();
    out << XmlWriter::attribute("start_atom") << object.start_atom();
    out << XmlWriter::attribute("connection") << object.connection();
  }


  void xml_handler< MrepConnection<2> >::read_object(XmlReader & in, MrepConnection<2> & object) const
  {
    float_t radius;
    float_t rotation;

    in >> XmlReader::attribute("r") >> radius;
    in >> XmlReader::attribute("rotation") >> rotation;

    object.assign(rotation, radius);
  }

  void xml_handler< MrepConnection<2> >::write_object(const MrepConnection<2> & object, XmlWriter & out) const
  {
    out << XmlWriter::attribute("r") << object.radius();
    out << XmlWriter::attribute("rotation") << object.rotation();
  }

  void xml_handler<MrepSkeleton2d>::read_object(XmlReader & in, MrepSkeleton2d & object) const
  {
    Position2d position;
    std::vector<MrepAtom> atoms;
    std::vector< MrepConnection<2> > connections;

    in >> position;
    in >> atoms;
    in >> connections;

    object.resize(atoms.size(), connections.size());
    object.set_position(position);
    
    for(std::size_t i = 0; i < atoms.size(); ++i)
      object.set_atom(i, atoms[i]);
    
    for(std::size_t i = 0; i < connections.size(); ++i)
      object.set_connection(i, connections[i]);
  }

  void xml_handler<MrepSkeleton2d>::write_object(const MrepSkeleton2d & object, XmlWriter & out) const
  {
    out << object.position();

    for(std::size_t i = 0; i < object.n_atoms(); ++i)
      out << object.atom(i);

    for(std::size_t i = 0; i < object.n_connections(); ++i)
      out << object.connection(i);
  }

  void xml_handler<MrepModel2d>::read_object(XmlReader & in, MrepModel2d & object) const
  {
    xml_handler<MrepSkeleton2d> skeleton_handler;
    
    skeleton_handler.read_object(in, object);
  }

  void xml_handler<MrepModel2d>::write_object(const MrepModel2d & object, XmlWriter & out) const
  {
    xml_handler<MrepSkeleton2d> skeleton_handler;
    
    skeleton_handler.write_object(object, out);
  }

  void xml_handler<PolygonModel2d>::read_object(XmlReader & in, PolygonModel2d & object) const
  {
    xml_handler<MrepSkeleton2d> skeleton_handler;
    
    skeleton_handler.read_object(in, object);
  }

  void xml_handler<PolygonModel2d>::write_object(const PolygonModel2d & object, XmlWriter & out) const
  {
    xml_handler<MrepSkeleton2d> skeleton_handler;
    
    skeleton_handler.write_object(object, out);
  }
}


