// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef XMLWRITER_H
#define XMLWRITER_H

#include <core/cio.hpp>

#include <boost/shared_ptr.hpp>


typedef struct _xmlDoc xmlDoc;
typedef xmlDoc * xmlDocPtr;

typedef struct _xmlTextWriter xmlTextWriter;
typedef xmlTextWriter * xmlTextWriterPtr;

namespace imaging
{
  template<class data_t>
  class xml_handler;
  
  /** \ingroup xml
      \brief Reads XML files using a stream interface.
      
      This class enables you to write an XML file using an interface which loosely resembles the stream interface used in \em iostream to write to output streams. An XmlWriter object is initialized by a root element and then represents an XML stream to which you can write using operator<<(const data_t &), where \em data_t is the type of the object to be written. XmlWriter provides the stream operator for the built-in data types. If called for other types, XmlWriter calls the function <em>xml_handler<data_t>::write_object()</em> which then writes the object to the stream. I.e. by specializing xml_handler for custom data types the user enables XmlWriter to write these types. Implementations of xml_handler for various classes of the \em imaging2 modules can be found in the files xmlio.hpp and xmlio.cxx in the module subdirectories.
      
      An XmlWriter object writes objects at the position of its <em>current path</em>. This path is the sequence of element names (and possibly indices) which leads to a given node in an XML file. This concept is very similiar to the idea of XPath. The current path of a XmlWriter can be altered by passing <em>stream commands</em> to the input stream. An example is given below:
  \code
  using imaging;
  
  XmlWriter xml_out("imaging2"); // construct an output stream; "imaging2" is the root element
  
  std::vector<float_t> some_floats; // a vector of floats
  
  xml_out << XmlWriter::element("some_float"); // set the current element to "some_float"
  xml_out << some_floats; // write the floats, each of them is enclosed by <some_float>...</some_float>
                          // there is no need to pass XmlWriter::end_element here,
                             because float is a built-in type
  
  std::string some_name("Franz");
  unsigned int some_age(50);
  
  xml_out << XmlWriter::element("person"); // add element "person" to the current path
  xml_out << XmlWriter::attribute("age") << some_age; // write attribute "age"
  xml_out << some_person; // write person's name
  xml_out << XmlWriter::end_element; // close element "person"
  
  xml_out.write_file("file.xml"); // write stream content to file
  \endcode
  Here, \em XmlWriter::element(), \em XmlWriter::attribute() and \em XmlWriter::end_element are stream commands.
  This code produces an XML file as the following:
  \code
  <imaging2>
    <person age="50">Franz</person>
    <some_float>...</some_float>
    <some_float>...</some_float>
  </imaging2>
  \endcode
  
  Assume further that the user has implemented a class \em Person and \em xml_handler<Person> which writes a person's data as above. Assuming that the class declaration and the declaration of \em xml_handler<Person> both reside in "Person.h", then she can use XmlWriter in the following way:
  \code
  #include "Person.h"
  
  XmlWriter xml_out("imaging2"); // construct an output stream; "imaging2" is the root element
  
  Person person;
  std::vector<Person> people;
  
  // fill person and people with data...
  
  xml_out << person; // write person (note that there is no need to add the
                     // element "person" to the current path)  
  xml_out << people; // write all the people at the current path
  \endcode
  This will result in one <tt>\<person/\></tt> element for the first person and one for each entry of the vector \em people.
      
      Note writing a vector of floats is different from a custom data type.
      For the custom data \em Person neither XmlWriter::element() nor XmlWriter::end_element have to be passed to the stream.
      In case of a vector of floats XmlWriter::element() must be passed, but not XmlWriter::end_element.
      The reason for this behavior is that XmlWriter does not know which element name identifies the floats. Introducing a standard element name (such as <tt>\<float\> \</float\></tt>) would often be inconvenient because it makes intuitive element names like <tt>\<super_parameter\>10^9\</super_parameter\></tt> impossible. This behavior is the same for all built-in types.
      
  \sa XmlReader
  */
  class XmlWriter
  {
    /** \cond */
    class stream_command
    {
      friend class XmlWriter;
      size_t _cmd_id;

    protected:
      enum commands { ELEMENT = 0, END_ELEMENT, ATTRIBUTE, COMMENT };

      std::string _identifier;

    public:
      stream_command() {}
      stream_command(size_t cmd_id, const std::string & identifier) : _cmd_id(cmd_id), _identifier(identifier) {}
    };

    template<size_t N, class data_t>
    void fixed_vector_2_string(const ublas::fixed_vector<data_t, N> & v, std::string & str)
    {
      std::ostringstream out;
      
      out << "(";
      for(std::size_t i = 0; i < N - 1; ++i)
        out << v(i) << ",";
        
      out << v(N - 1) << ")";
        
      str = out.str();
    }
    
    template <class data_t>
    XmlWriter &  write_integral_type_value(const data_t & value);
    
    template <class data_t>
    XmlWriter &  write_integral_type_vector(const std::vector<data_t> & vector);
    
    /** \endcond */
    
  public:
    /** \brief Stream command to add an %element to the current path. */ 
    class element : public stream_command
    {
    public:
      /** Passing this element to an XML output stream adds the element \em identifier to the current path. */
      element(const std::string & identifier) : stream_command(ELEMENT, identifier) {}
    };
    /** \brief Stream command to add an %attribute at the current path. */ 
    class attribute : public stream_command
    {
    public:
      attribute(const std::string & identifier) : stream_command(ATTRIBUTE, identifier) {}
    };
    
    /** \brief Stream command to add a %comment at the current path. */
    class comment : public stream_command
    {
    public:
        /** Passing this element to an XML output stream adds the comment \em identifier at the current path. */
      comment(const std::string & identifier) : stream_command(COMMENT, identifier) {}
    };

  private:
    xmlDocPtr _doc;
    xmlTextWriterPtr _writer;

    bool _attribute_active;
    std::string _current_attribute;
    std::string _current_attribute_content;
    
    bool _element_active;
    std::string _current_element;

    void flush_element();

    static const char ISO_8859_1 [];


  public:

    /** \brief Stream command to end an %element. */
    static const stream_command end_element;

    /** Constructs an XmlWriter object and adds \em root_element. If \em style_file_name contains the address of an XSL file, an \c xml-stylesheet tag pointing to the the file is added to the resulting XML file. */ 
    XmlWriter(const std::string & root_element, const std::string & style_file_name = "");

    ~XmlWriter();

    /** Passes any other stream command to the stream. */
    XmlWriter & operator<<(const stream_command & command);
    
    /** Passes an XmlReader::element() command to the stream. */
    XmlWriter & operator<<(const element & command);
    
    /** Passes an XmlReader::attribute() command to the stream. */
    XmlWriter & operator<<(const attribute & command);
    
    /** Writes a comment to the stream. */
    XmlWriter & operator<<(const comment & command);
    
    /** Writes a string to the stream. */
    XmlWriter & operator<<(const std::string & string);
    
    /** Writes a character to the stream. */
    XmlWriter & operator<<(const char* str);
    
    /** Writes a floating point value to the stream. */
    XmlWriter & operator<<(float_t value);
    
    /** Writes a unsigned integer  to the stream. */
    XmlWriter & operator<<(size_t value);

    /** Writes a vector of tla::vector objects to the stream. */    
    template<size_t N, class data_t>
    XmlWriter & operator<<(const ublas::fixed_vector<data_t, N> & value)
    {
      std::string temp;
      fixed_vector_2_string(value, temp);
      
      *this << temp;
    
      return *this;
    }
    
    /** Writes a \em data_t object to the stream. The class specialization xml_handler<data_t>::write_object() has to be defined. */
    template<class data_t>
    XmlWriter & operator<<(const data_t & object)
    {
      xml_handler<data_t> handler;
    
      *this << element(handler.element_name);
      handler.write_object(object, *this);
      *this << end_element;

      return *this;
    }
    
    /** Writes a vector of \em data_t objects to the stream. The class specialization xml_handler<data_t>::write_object() has to be defined. */   
    template<class data_t>
    XmlWriter & operator<<(const std::vector<data_t> & vector)
    {
      for(size_t i = 0; i < vector.size(); ++i)
        *this << vector[i];

      return *this;
    }
    
    /** Writes a vector of floats objects to the stream. Each float is enclosed by element tags defined by the \em last XmlWriter::element() passed to the input stream. In this respect this function is different from the generic version operator<<(data_t &). */   
    XmlWriter & operator<<(const std::vector<float_t> & vector);
    
    /** Writes a vector of \em data_t objects to the stream. The class specialization xml_handler<data_t>::write_object() has to be defined. */   
    template<class data_t>
    XmlWriter & operator<<(const std::vector< boost::shared_ptr<data_t> > & vector)
    {
      for(size_t i = 0; i < vector.size(); ++i)
        *this << *vector[i];

      return *this;
    }

    /** Writes the content of the stream to \em file_name. */
    void write_file(const std::string & file_name);
  };
}

#endif
