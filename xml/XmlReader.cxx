#include <xml/XmlReader.hpp>

#include <core/MessageInterface.hpp>

#include <iostream>
#include <sstream>

#include <libxml/xpathInternals.h>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>



namespace imaging
{
  const XmlReader::stream_command XmlReader::end_element(stream_command::END_ELEMENT, "");
  
  size_t XmlReader::n_elements() const
  {
    std::string path;
  
    compute_current_path(path);
    return compute_n_elements(path);
  }
  
  size_t XmlReader::n_child_elements() const
  {
    std::string path;
    
    compute_current_path(path);
    xmlNodeSetPtr node_set_ptr;
    size_t n_nodes = evaluate_x_path(path, node_set_ptr);
    
    size_t n = 0;
    for(size_t i = 0; i < n_nodes; ++i)
    {
      xmlNode * child_node = &(node_set_ptr->nodeTab[i]->children[0]);
      while(child_node->next != 0x0)
      {
        if(child_node->type == XML_ELEMENT_NODE)
          ++n;
        
        child_node = child_node->next;
      }
    }
    
    return n;  
  }

  XmlReader & XmlReader::operator>>(const stream_command & command)
  {
    std::string element_name;
    
    switch(command._cmd_id)
    {
    case stream_command::ATTRIBUTE:
      _current_attribute = command._identifier;
      _attribute_active = true;
      return *this;
      break;
      
    case stream_command::ELEMENT:
    case stream_command::INDEX:
    case stream_command::ELEMENT_WITH_INDEX:
    case stream_command::CHILD_ELEMENT:
      _current_path.push_back(command);
      _attribute_active = false;
      return *this;
      break;
      
    case stream_command::END_ELEMENT:
      _current_path.pop_back();
      break;
    }

    return *this;
  }

  XmlReader & XmlReader::operator>>(std::string & string)
  {
    if(! _attribute_active)
    {
      std::string path;
      compute_current_path(path);

      parse_element(path, string);
    }
    else
    {
      _attribute_active = false;
      
      parse_attribute(_current_attribute, string);
    }

    return *this;
  }
  

  XmlReader::XmlReader(const std::string & xml_file_name, const std::string & root_element) : _current_attribute(""), _attribute_active(false)
  {
    xmlInitParser();

    _doc = xmlParseFile(xml_file_name.c_str());
    if (_doc == NULL)
      throw FileIoException("FileIoException: unable to parse file " + xml_file_name + ".");

    _xpathCtx = xmlXPathNewContext(_doc);

    if(_xpathCtx == 0)
      throw Exception("Exception: unable to create new XPath context.");

    MessageInterface::out("Read parameter file " + xml_file_name + ".", MessageInterface::DEBUG_ONLY);

    _current_attribute = "";
    _attribute_active = false;

    *this >> element(root_element);
  }


  XmlReader::~XmlReader()
  {
    if(_xpathCtx)
      xmlXPathFreeContext(_xpathCtx);
    if(_doc)
      xmlFreeDoc(_doc);
  }


  size_t XmlReader::evaluate_x_path(const std::string & path_expr, xmlNodeSetPtr & node_set_ptr) const
  {
    xmlXPathObjectPtr xpathObj;

    xpathObj = xmlXPathEvalExpression((xmlChar*)path_expr.c_str(), _xpathCtx);
    if(xpathObj == 0)
      throw Exception("Could not evaluate path '" + path_expr + "'.");

    node_set_ptr = xpathObj->nodesetval;

    if(node_set_ptr == 0x0)
      throw Exception("Could not evaluate nodeset of path '" + path_expr + "'.");

    return node_set_ptr->nodeNr;
  }


  void XmlReader::parse_element(const std::string & path_expr, std::string & result) const
    throw (XmlNoTagException)
  {
    xmlNodeSetPtr node_set_ptr;
    size_t n_elements;

    n_elements = evaluate_x_path(path_expr, node_set_ptr);

    if(n_elements == 0)
      throw XmlNoTagException("No element '" + path_expr + "' in XmlReader::parse_element().");
    else
    {
      xmlNodePtr node_ptr = node_set_ptr->nodeTab[0];

      char * content = (char*)xmlNodeListGetString(_doc, node_ptr->children, 1);
      if (content == 0x0)
      {
        //  MessageInterface::out("No content in element '" + path_expr + "'.", MessageInterface::DEBUG, 1);
        //  throw XmlNoTagException();
        result = "";
      }
      else
        result = content;

      MessageInterface::out("Read element '" + path_expr + "': " + result + ".", MessageInterface::DEBUG_ONLY, 1);

      xmlFree(content);
    }
  }


  void XmlReader::parse_attribute(const std::string & attribute, std::string & result) const
    throw (XmlNoTagException)
  {
    std::string path_expr;
    size_t n_elements;
    xmlNodeSetPtr node_set_ptr;
    xmlNodePtr node_ptr;

    compute_current_path(path_expr);
    n_elements = evaluate_x_path(path_expr, node_set_ptr);

    if(n_elements == 0)
      throw XmlNoTagException("No element '" + path_expr + "' (while parsing attribute '" + attribute + "') in XmlReader::parse_attribut().");

    node_ptr = node_set_ptr->nodeTab[0];

    char * content = (char*)xmlGetProp(node_ptr, (xmlChar*)attribute.c_str());
    if (content == 0x0)
      throw XmlNoTagException("No attribute '" + attribute + "' in element '" + path_expr + "' in XmlReader::parse_attribute().");
    else
    {
      result = content;
      MessageInterface::out("Read attribute '" + attribute + "' of element '" + path_expr + "': " + result + ".", MessageInterface::DEBUG_ONLY, 1);
    }

  }

  void XmlReader::compute_current_path(std::string & path) const
  {
    path = "";
    
    for(size_t i = 0; i < _current_path.size(); ++i)
    {
      switch(_current_path[i]._cmd_id)
      {
      case stream_command::ELEMENT:
        path += "/" + _current_path[i]._identifier;
        break;
        
      case stream_command::ELEMENT_WITH_INDEX:
        path += "/" + _current_path[i]._identifier + "[" + boost::lexical_cast<std::string>(_current_path[i]._index + 1) + "]";
        break;
        
      case stream_command::CHILD_ELEMENT:
        path += "/*[" + boost::lexical_cast<std::string>(_current_path[i]._index + 1) + "]";
        break;
        
      case stream_command::INDEX:
        path += "[" + boost::lexical_cast<std::string>(_current_path[i]._index + 1) + "]";
        break;
      }
    }     
  }

  size_t XmlReader::compute_n_elements(const std::string & path_expr) const
  {
    xmlNodeSetPtr node_set_ptr;

    return evaluate_x_path(path_expr, node_set_ptr);
  }
  
  std::string XmlReader::current_element_name() const
  {
    xmlNodeSetPtr node_set_ptr;
    std::string path;
    size_t n_elements;
    
    compute_current_path(path);
    n_elements = evaluate_x_path(path, node_set_ptr);
    
    if(n_elements == 0)
      throw XmlNoTagException("No element '" + path + "' in XmlReader::parse_element().");
    
    xmlNodePtr node_ptr = node_set_ptr->nodeTab[0];
    
    return std::string((char*)node_ptr->name);
  }
  
  XmlReader & XmlReader::operator>>(std::vector<float_t> & vector)
  {
    return read_integral_type_vector(vector);
  }
  
  XmlReader & XmlReader::operator>>(float_t & value)
  {
    return read_integral_type_value(value);
  }
  
  XmlReader & XmlReader::operator>>(size_t & value)
  {
    return read_integral_type_value(value);
  }
  
  XmlReader & XmlReader::operator>>(int & value)
  {
    return read_integral_type_value(value);
  }
  
  XmlReader & XmlReader::operator>>(bool & value)
  {
    return read_integral_type_value(value);
  }
  
  XmlReader & XmlReader::operator>>(std::string & string);
  
  XmlReader & XmlReader::operator>>(std::vector<std::string> & vector)
  {
    return read_integral_type_vector<std::string>(vector);
  }


}

