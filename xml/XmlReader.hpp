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

#ifndef XML_XMLREADER_H
#define XML_XMLREADER_H

#include <core/cio.hpp>
#include <core/Exception.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

typedef struct _xmlDoc xmlDoc;
typedef xmlDoc * xmlDocPtr;

typedef struct _xmlXPathContext xmlXPathContext;
typedef xmlXPathContext * xmlXPathContextPtr;

typedef struct _xmlNodeSet xmlNodeSet;
typedef xmlNodeSet * xmlNodeSetPtr;
  
namespace imaging
{
  
  /** \ingroup xml
      \brief Utility class template for XmlReader and XmlWriter.
      
      XmlReader and XmlWriter will call 
  \code
  void xml_handler<data_t>::read_object(XmlReader & in, data_t & object) const;
  void xml_handler<data_t>::write_object(const data_t & object, XmlWriter & out) const;
  \endcode
      if they encounter a \em data_t object as argument of the input or output stream operator. Thus, to support XML input and/or output for a custom data type, the user has to specialize xml_handler and implement the above functions.
  */
  template<class data_t>
  class xml_handler;
  
  /** \ingroup xml
      \brief Reads XML files using a stream interface.
      
      This class enables you to read an XML file using an interface which loosely resembles the stream interface used in \em iostream to read from input streams. An XmlReader object is initialized by an XML file and then represents an XML stream which you can read using operator>>(data_t &), where \em data_t is the type of the object to be read. XmlReader provides the stream operator for the built-in data types. If called for other types, XmlReader calls the function <em>xml_handler<data_t>::read_object()</em> which then reads the object from the stream.  I.e. by specializing xml_handler for custom data types the user enables XmlReader to read these types. Implementations of xml_handler for various classes of the \em imaging2 modules can be found in the files xmlio.hpp and xmlio.cxx in the module subdirectories.
       
      Consider the following XML file:
  \code
  <imaging2>
    <person age="50">Franz</person>
    <person age="50">Arne</person>
    
    <some_float>0.5</some_float>
    <some_float>1.0</some_float>
    <some_float>1.5</some_float>
  </imaging2>
  \endcode
      
      An XmlReader object reads objects at the position of its <em>current path</em>. This path is the sequence of element names (and possibly indices) which leads to a given node in an XML file. This concept is very similiar to the idea of XPath. The current path of a XmlReader can be altered by passing <em>stream commands</em> to the input stream. An example is given below:
  \code
  using imaging;
  
  XmlReader xml_in("some_file.xml", "imaging2"); // read some_file.xml; "imaging2" is the root element
  
  std::string name;
  unsigned int age = 0;
  
  xml_in >> XmlReader::element("person"); // select element "person"
  xml_in >> XmlReader::attribute("age") >> age; // read attribute "age"
  xml_in >> name >> XmlReader::end_element; // read the person's name and close element "person"
  
  std::vector<float_t> some_floats;
  xml_in >> XmlReader::element("some_float"); // select element "some_float"
  xml_in >> some_floats; // read the floats                         
                         // there is no need to pass XmlReader::end_element here,
                            because float is a built-in type
  
  \endcode
     Here, \em XmlReader::element(), \em XmlReader::attribute() and \em XmlReader::end_element are stream commands.
      In case of the above XML file, the vector \em some_floats will be resized to length 3 and will hold the data <em>0.5, 1.0, 1.5</em>.
  
      Assume further that the user has implemented a class \em Person and \em xml_handler<Person> which reads a person's data as above. Assuming that the class declaration and the declaration of \em xml_handler<Person> both reside in "Person.h", then she can use XmlReader in the following way:
  \code
  #include "Person.h"
  
  XmlReader xml_in("some_file.xml", "imaging2"); // read some_file.xml; "imaging2" is the root element
  
  Person person;
  std::vector<Person> people;
  
  xml_in >> person; // read person (note that there is no need to select the element "person")
  xml_in >> people; // read all the people at the current path
  \endcode
      In case of the above XML file, the vector \em people will be resized to length 1 and will hold the data of fifty years old Franz.
      
      Note reading a vector of floats is different from a custom data type.
      For the custom data \em Person neither XmlReader::element() nor XmlReader::end_element have to be passed to the stream.
      In case of a vector of floats XmlReader::element() must be passed, but not XmlReader::end_element.
      The reason for this behavior is that XmlReader does not know which element name identifies the floats. Introducing a standard element name (such as <tt>\<float\> \</float\></tt>) would often be inconvenient because it makes intuitive element names like <tt>\<super_parameter\>10^9\</super_parameter\></tt> impossible. This behavior is the same for all built-in types.
      
  \sa XmlWriter
  */
  class XmlReader
  {
    /** \cond */
    class stream_command
    {
      friend class XmlReader;

    protected:
      enum commands { ELEMENT = 0, INDEX, ELEMENT_WITH_INDEX, END_ELEMENT, ATTRIBUTE, CHILD_ELEMENT };

      size_t _cmd_id;
      std::string _identifier;
      size_t _index;

    public:
      stream_command() {}
      stream_command(size_t cmd_id, const std::string & identifier, size_t index = 0) : 
        _cmd_id(cmd_id), _identifier(identifier), _index(index) {}}
    ;
    
  /** \endcond */

  public:
    /** \brief Stream command to add an %element to the current path. */ 
    class element : public stream_command
    {
    public:
      /** Passing this element to an XML input stream adds the element \em identifier to the current path. */ 
      element(const std::string & identifier) : stream_command(ELEMENT, identifier) {}
      
      /** Passing this element to an XML input stream adds the \em index-th (zero-based) element named \em identifier to the current path. */ 
      element(const std::string & identifier, size_t index) :
        stream_command(ELEMENT_WITH_INDEX, identifier, index) {}
      
      /** Passing this element to an XML input stream selects the \em index-th (zero-based) element of the elements selected by the current path. */ 
      element(size_t index) : stream_command(INDEX, "", index) {}
    };
    
    /** \brief Stream command to add a child %element to the current path. */ 
    class child_element : public stream_command
    {
    public:
      /** Passing this element to an XML input stream adds \em index-th child %element of the current %element to the current path. */ 
      child_element(size_t index) : stream_command(CHILD_ELEMENT, "", index) {}
    };

    
    /** \brief Stream command to select an %attribute at the current path. */  
    class attribute : public stream_command
    {
    public:
      /** Passing this attribute to an XML input stream selects the %attribute \em identifier at the current path. */ 
      attribute(const std::string & identifier) : stream_command(ATTRIBUTE, identifier) {}}
    ;

    /** \brief Stream command to read an object if it is present at the current path and to retrieve an default value if not.
    
        Assume that the user wants to read an object from an XML input stream but she is not sure if it can actually be found in the stream. Assume further that in case the object can \em not be read the user just wants to set it to a default value. Then she can use a default_value to accomplish this task in a compact way.
    
        The alternate way to handle this situation would be to try to read the object and to explicitly handle XmlNoTagException exceptions.
    */
    template <class data_t>
    class default_value
    {
      friend class XmlReader;

      data_t & _object;
      const data_t & _default_object;

    public:
      /** Passing this command to an XML input stream reads \em object from the stream if it is present at the current stream path. If not, it sets \em object to \em default_object.
      */
      default_value(data_t & object, const data_t & default_object) : _object(object), _default_object(default_object) {}
    };
    
    /** \brief Stream command to read an object if it is present at the current path and to keep it unchanged if not.
    
        Assume that the user wants to read an object from an XML input stream but she is not sure if it can actually be found in the stream. Assume further that in case the object can \em not be read the user just wants the object to be the same as before without any further error handling. Then she can use a try_read to accomplish this task in a compact way.
    
        The alternate way to handle this situation would be to try to read the object and to explicitly handle XmlNoTagException exceptions.
    */
    template <class data_t>
    class try_read
    {
      friend class XmlReader;

      data_t & _object;

    public:
      /** Passing this command to an XML input stream reads \em object from the stream if it is present at the current stream path. If not, it keeps \em object unchanged.
      */
      try_read(data_t & object) : _object(object){}
    };


  /** \brief This exception is throw if the attempt to read an object fails. 
  
  An XmlNoTagException is thrown whenever the user attempts to read an object, which is not present at the current path. This does not hold when reading a vector of objects. In this case the length of the vector is simply set to 0, if no object can be read. Also note that an XmlNoTagException exception is \em never thrown if the user passes a stream command to the input stream. */
  class XmlNoTagException : public Exception
    {
    public:
      /** Constructor. The user can pass a custom error message. */
      XmlNoTagException(std::string msg = std::string("XmlNoTagException!")) : Exception(msg) {}}
    ;

  private:
    xmlDocPtr _doc;
    xmlXPathContextPtr _xpathCtx;

    std::vector<stream_command> _current_path;
    std::string _current_attribute;
    bool _attribute_active;

    void parse_element(const std::string & path_expression, std::string & result) const
      throw (XmlNoTagException);
    
    void parse_attribute(const std::string & attribute, std::string & result) const
      throw (XmlNoTagException);

    size_t evaluate_x_path(const std::string & path_expr, xmlNodeSetPtr & node_set_ptr) const;

    void compute_current_path(std::string & path) const;

    size_t compute_n_elements(const std::string & path_expr) const;

    template<size_t N, class data_t>
    void string_2_fixed_vector(const std::string & str, ublas::fixed_vector<data_t, N> & v)
    {
      std::istringstream in(str);
      char buffer[5];

      in.getline(buffer, 5, '(');

      for(std::size_t i = 0; i < N - 1; ++i)
      {
        v(i) = data_t(float_t(0));
        in >> v(i);
        in.getline(buffer, 5, ',');
      }

      v(N - 1) = data_t(float_t(0));
      in >> v(N - 1);
    }
    
    template <class data_t>
    XmlReader &  read_integral_type_value(data_t & value)
    {
      std::string temp;
  
      *this >> temp;
  
      try
      {
        value = boost::lexical_cast<data_t>(temp);
      }
      catch(boost::bad_lexical_cast & e)
      {
        throw Exception("Exception: lexical cast of '" + temp + "' failed in XmlReader::read_integral_type_value().");
      }
  
      return *this;
    }
    
    template <class data_t>
    XmlReader & read_integral_type_vector(std::vector<data_t> & vector)
    {
      std::string path;
      size_t n_elements;
  
      compute_current_path(path);
      n_elements = compute_n_elements(path);
      vector.resize(n_elements);
  
      for(size_t i = 0; i < n_elements; ++i)
      {
        *this >> element(i);
        *this >> vector[i];
        *this >> end_element;
      }
      
      *this >> end_element;
  
      return *this;
    }

  public:
    /** \brief Stream command to end an element.
    
        Stream command to remove the last element from the current path. Every XmlReader::element(std::string &) and XmlReader::element(std::string &, size_t) command passed to the input stream has to be followed by an Xml::end_element command. This is not true for XmlReader::element(size_t) and XmlReader::attribute() commands. */
    static const stream_command end_element;

    /** Constructs an XML input stream from file \em file_name and sets the current path to \em root_element. An XML file can only have one root node, so it is okay to fix it during construction. */
    XmlReader(const std::string & xml_file_name, const std::string & root_element);

    /** Destructor. Closes the input file. */
    ~XmlReader();

    /** Passes any other stream command to the stream. */
    XmlReader & operator>>(const stream_command & command);

    /** Reads a ublas::fixed_vector from the stream. */
    template<size_t N, class data_t>
    XmlReader & operator>>(ublas::fixed_vector<data_t, N> & value)
    {
      std::string temp;

      (*this) >> temp;

      string_2_fixed_vector(temp, value);

      return *this;
    }

    /** Reads a vector of fixed size vectors from the stream. The function assumes that each vector is enclosed by element tags identified by the \em last XmlReader::element() passed to the input stream. In this respect this function is different from the generic version operator>>(data_t &). */
    template<size_t N, class data_t>
    XmlReader & operator>>(std::vector< ublas::fixed_vector<data_t, N> > & vector)
    {
      std::string path;
      size_t n_elements;

      compute_current_path(path);
      n_elements = compute_n_elements(path);
      vector.resize(n_elements);

      for(size_t i = 0; i < n_elements; ++i)
      {
        *this >> element(i);
        *this >> vector[i];
        *this >> end_element;
      }
      
      *this >> end_element;

      return *this;
    }
    
    /** Reads a vector of floats from the stream. The function assumes that each string is enclosed by element tags identified by the \em last XmlReader::element() passed to the input stream. In this respect this function is different from the generic version operator>>(data_t &). */
    XmlReader & operator>>(std::vector<float_t> & vector);
    
    /** Reads a floating point value from the stream. */
    XmlReader & operator>>(float_t & value);
    
    /** Reads an unsigned integer from the stream. */
    XmlReader & operator>>(size_t & value);
    
    /** Reads an integer from the stream. */
    XmlReader & operator>>(int & value);
    
    /** Reads a bool from the stream. */
    XmlReader & operator>>(bool & value);
    
    /** Reads string from the stream. */
    XmlReader & operator>>(std::string & string);
    
    /** Reads a vector of strings from the stream. The function assumes that each string is enclosed by element tags identified by the \em last XmlReader::element() passed to the input stream. In this respect this function is different from the generic version operator>>(data_t &). */
    XmlReader & operator>>(std::vector<std::string> & vector);

    /** Passes an XmlReader::default_value() command to the stream. */
    template<class data_t>
    XmlReader & operator>>(const default_value<data_t> & command)
    {
      try { *this >> command._object; }
      catch (XmlNoTagException & exception) { command._object = command._default_object; }

      return *this;
    }

    /** Passes an XmlReader::try_read() command to the stream. */
    template<class data_t>
    XmlReader & operator>>(const try_read<data_t> & command)
    {
      try { *this >> command._object; }
      catch (XmlNoTagException & exception) { }

      return *this;
    }
    
    /** Reads a vector of \em data_t objects from the stream. The class specialization xml_handler<data_t> has to be defined. */
    template<class data_t>
    XmlReader & operator>>(std::vector<data_t> & vector)
    {
      xml_handler<data_t> handler;
      
      std::string path;
      size_t n_elements;
      
      *this >> element(handler.element_name);

      compute_current_path(path);
      n_elements = compute_n_elements(path);
      vector.resize(n_elements);

      for(size_t i = 0; i < n_elements; ++i)
      {
        *this >> element(i);
        handler.read_object(*this, vector[i]);
        *this >> end_element;
      }

      *this >> end_element;

      return *this;
    }
    
    /** Constructs \em data_t objects in \em vector and reads \em data_t objects from the stream. The class specialization xml_handler<data_t> has to be defined. */
    template<class data_t>
    XmlReader & operator>>(std::vector< boost::shared_ptr<data_t> > & vector)
    {
      xml_handler<data_t> handler;
      
      std::string path;
      size_t n_elements;
      
      *this >> element(handler.element_name);

      compute_current_path(path);
      n_elements = compute_n_elements(path);
      vector.resize(n_elements);
      for(size_t i = 0; i < n_elements; ++i)
      {
        try
        {
          vector[i].reset(new data_t);
        }
        catch(std::bad_alloc & e)
        {
          n_elements = 0;
          vector.resize(n_elements);
        }
      }

      for(size_t i = 0; i < n_elements; ++i)
      {
        *this >> element(i);
        handler.read_object(*this, *vector[i]);
        *this >> end_element;
      }

      *this >> end_element;

      return *this;
    }
    
    /** Reads a \em data_t object from the stream. The class specialization xml_handler<data_t>::read_object() has to be defined. */
    template<class data_t>
    XmlReader & operator>>(data_t & object)
    {
      xml_handler<data_t> handler;
      bool xml_exception_thrown = false;
      XmlReader::XmlNoTagException exception;
      
      *this >> element(handler.element_name);
      
      try
      {
        handler.read_object(*this, object);
      }
      catch(XmlReader::XmlNoTagException & e)
      {
        xml_exception_thrown = true;
        exception = e;
      }
      
      *this >> end_element;
      
      if(xml_exception_thrown)
        throw exception;

      return *this;
    }

    /** Returns the number of different nodes in the document which are referenced by the current path. If the current path uniquely determines an object this function returns 1. If no object exists at the current path it returns 0. */
    size_t n_elements() const;

    /** Returns the number of child nodes of elements in the document which are referenced by the current path. */
    size_t n_child_elements() const;

    /** Returns the name of the current element. In case the last element was selected by passing <em>XmlReader::element(identifier)</em> to the stream, this function returns \em identifier. */
    std::string current_element_name() const;
  };
}

#endif
