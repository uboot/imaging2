
#include <core/Cmessage.hpp>
#include <xml/XmlReader.hpp>
#include <xml/XmlWriter.hpp>

using namespace imaging;

class Person
{
  std::string _name;
  std::size_t _age;
  
public:
  const std::string & name() const { return _name; }
  std::size_t age() const { return _age; }
  
  void set_name(const std::string & name) { _name = name; }
  void set_age(const std::size_t & age) { _age = age; }
};


namespace imaging
{
  template<>
  class xml_handler<Person>
  {
  public:
    static const std::string element_name;
    
    void read_object(XmlReader & in, Person & object) const
    {
      std::string name;
      std::size_t age;
      
      in >> XmlReader::attribute("age") >> XmlReader::default_value<std::size_t>(age, 0);
      in >> name;
      
      object.set_name(name);
      object.set_age(age);
    }
    
    void write_object(const Person & object, XmlWriter & out) const
    {
      out << XmlWriter::attribute("age") << object.age();
      out << object.name();
    }
  };
  
  const std::string xml_handler<Person>::element_name = "person";
}

int main ( int argc, char **argv )
{
  try
  {
    img::size_t verbose_level = Cmessage::ALL_MESSAGES;
    img::size_t main_argument_base_index = 1;
    
    if(argc < 2)
      throw Exception("Too little command line parameters provided!");
      
    if(argc == 3)
    {
      verbose_level = boost::lexical_cast<img::size_t>(argv[1]);
      main_argument_base_index = 2;
    }
    
    Cmessage::out.set_verbose_level(verbose_level);
    
    XmlReader xml_in(argv[main_argument_base_index], "imaging");
    std::vector< boost::shared_ptr<Person> > people;
    std::vector<img::float_t> floats;
    
    xml_in >> XmlReader::element("some_float");
    xml_in >> floats;
    xml_in >> people;
    
    for(std::vector<img::float_t>::iterator iter = floats.begin(); iter != floats.end(); ++iter)
      std::cout << *iter << std::endl;
      
    for(std::vector< boost::shared_ptr<Person> >::iterator iter = people.begin(); iter != people.end(); ++iter)
      (*iter)->set_age((*iter)->age() + 10);
    
    XmlWriter xml_out("imaging");
    
    xml_out << XmlWriter::element("some_float");
    xml_out << floats;
    xml_out << people;
    
    xml_out.write_file("out.xml");
  }

  catch ( Exception &exception )
  {
    std::cerr << exception.error_msg() << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
