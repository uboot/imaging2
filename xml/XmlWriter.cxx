#include <xml/XmlWriter.hpp>

#include <core/MessageInterface.hpp>

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include <boost/lexical_cast.hpp>


namespace imaging
{
  template <class data_t>
  XmlWriter &  XmlWriter::write_integral_type_value(const data_t & value)
  {
    *this << boost::lexical_cast<std::string>(value);
  
    return *this;
  }
  
  template <class data_t>
  XmlWriter &  XmlWriter::write_integral_type_vector(const std::vector<data_t> & vector)
  {
    std::string element_name = _current_element;
    
    if(vector.size() == 0)
    {
      _element_active = false;
      return *this;
    }
    
    typename std::vector<data_t>::const_iterator iter = vector.begin();
    *this << *iter << end_element;
    ++iter;
    
    for(; iter != vector.end(); ++iter)
      *this << element(_current_element) << *iter << end_element;
  
    return *this;
  }
  
  const XmlWriter::stream_command XmlWriter::end_element(stream_command::END_ELEMENT, "");

  const char XmlWriter::ISO_8859_1 [] = "ISO-8859-1";

  XmlWriter::XmlWriter(const std::string & root_element, const std::string & style_file_name)
  {
    int rc;

    /* Create a new XmlWriter for DOM, with no compression. */
    _writer = xmlNewTextWriterDoc(&_doc, 0);
    if (_writer == NULL)
      throw Exception("Exception: Error creating the xmlwriter.");


    /* Start the document with the xml default for the version,
     * encoding ISO 8859-1 and the default for the standalone
     * declaration. */
    rc = xmlTextWriterStartDocument(_writer, NULL, ISO_8859_1, NULL);
    if (rc < 0)
      throw Exception("Exception: Error at xmlTextWriterStartDocument.");

    if(style_file_name != "")
    {
      std::string style_file_tag = "type=\"text/xsl\" href=\"" + style_file_name + "\"";
      rc = xmlTextWriterWriteProcessingInstruction
           (_writer, BAD_CAST "xml-stylesheet", BAD_CAST style_file_tag.c_str());
      if (rc < 0)
        throw Exception("Exception: Error at xmlTextWriterWriteProcessingInstruction.");
    }

    _attribute_active = false;
    _current_attribute = "";
    _element_active = false;
    *this << element(root_element);
  }

  XmlWriter::~XmlWriter()
  {
    if(_writer)
      xmlFreeTextWriter(_writer);
  }

  XmlWriter & XmlWriter::operator<<(const element & command)
  {
    flush_element();
    
    _element_active = true;
    _current_element = command._identifier;
    
    return *this;
  }

  XmlWriter & XmlWriter::operator<<(const attribute & command)
  {
    flush_element();

    _attribute_active = true;
    _current_attribute = command._identifier;

    return *this;
  }

  void XmlWriter::flush_element()
  {
    if(_element_active)
    {
      int rc;
  
      rc = xmlTextWriterStartElement(_writer, BAD_CAST _current_element.c_str());
      if (rc < 0)
        throw Exception("Exception: Error at xmlTextWriterStartElement.");
      _element_active = false;
    }
    
    if(_attribute_active)
    {
      int rc;

      rc = xmlTextWriterWriteAttribute(_writer, BAD_CAST _current_attribute.c_str(),
                                       BAD_CAST _current_attribute_content.c_str());
      if (rc < 0)
        throw Exception("Exception: Error at xmlTextWriterWriteAttribute.");

      _attribute_active = false;
      _current_attribute = "";
      _current_attribute_content = "";
    }
  }



  XmlWriter & XmlWriter::operator<<(const stream_command & command)
  {
    flush_element();

    int rc;

    switch(command._cmd_id)
    {
    case stream_command::END_ELEMENT:
      rc = xmlTextWriterEndElement(_writer);
      if (rc < 0)
        throw Exception("Exception: Error at xmlTextWriterEndElement.");
      *this << "\n";
      break;
    }

    return *this;
  }
  
  XmlWriter & XmlWriter::operator<<(const comment & command)
  {
    flush_element();

    int rc;

    rc = xmlTextWriterWriteComment(_writer, BAD_CAST command._identifier.c_str());
    if (rc < 0)
      throw Exception("Exception: Error at xmlTextWriterWriteComment.");

    return *this;

  }

  XmlWriter & XmlWriter::operator<<( const std::string & string )
  {
    if(_attribute_active)
    {
      _current_attribute_content += string;
      flush_element();
    }
    else
    {
      flush_element();
      
      int rc;

      rc = xmlTextWriterWriteFormatString(_writer, "%s",
                                          BAD_CAST string.c_str());
      if (rc < 0)
        throw Exception("Exception: Error at xmlTextWriterWriteFormatString.");
    }

    return *this;
  }

  void XmlWriter::write_file( const std::string & output_file_name )
  {
    *this << end_element; 
    
    int rc;
    rc = xmlTextWriterEndDocument(_writer);
    if (rc < 0)
      throw Exception("Exception: Error at xmlTextWriterEndDocument.");

    if(_writer)
      xmlFreeTextWriter(_writer);
    _writer = 0;

    rc = xmlSaveFileEnc(output_file_name.c_str(), _doc, ISO_8859_1);
    if (rc < 0)
      throw FileIoException("FileIoException: Error at xmlTextWriterEndDocument.");

    MessageInterface::out("Wrote file " + output_file_name + ".", MessageInterface::DEBUG_ONLY);
  }
  
  XmlWriter & XmlWriter::operator<<(const char* str) 
  {
    return write_integral_type_value(str);
  }
  
  XmlWriter & XmlWriter::operator<<(float_t value)
  {
    return write_integral_type_value(value);
  }
  
  XmlWriter & XmlWriter::operator<<(size_t value)
  {
    return write_integral_type_value(value);
  }
 
  XmlWriter & XmlWriter::operator<<(const std::vector<float_t> & vector)
  { 
    return write_integral_type_vector(vector);
  }
}

