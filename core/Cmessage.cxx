#include <core/Cmessage.hpp>

#include <core/utilities.hpp>

namespace imaging
{
  Cmessage Cmessage::out;
  MessageInterface & MessageInterface::out(Cmessage::out);
  
  Cmessage::Cmessage()
  {
    _current_intendation_level = 0; 
    _current_verbose_level = ALL_MESSAGES;
  }

  void Cmessage::set_verbose_level(size_t level)
  {
    _current_verbose_level = level;
  }

  void Cmessage::operator()(const std::string & message, const int priority_level, int intend)
  {
    (*this)(intend);

    if(priority_level >= _current_verbose_level)
    {

      for(int i = 0; i < _current_intendation_level; ++i)
        std::cout << "  ";

      std::cout << message << std::endl;
    }
    
    (*this)(-intend);
  }

  void Cmessage::operator()(int intend)
  {
    _current_intendation_level += intend;
    _current_intendation_level = max(_current_intendation_level, 0);
  }
}

