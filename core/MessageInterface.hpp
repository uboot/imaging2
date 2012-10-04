// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


#ifndef CORE_MESSAGEINTERFACE_H
#define CORE_MESSAGEINTERFACE_H

#include <core/imaging2.hpp>
#include <iostream>

namespace imaging
{

  /** \ingroup core
      \brief Abstract class interface for displaying log messages.
       
      This class defines the \em imaging2 log message interface. The user can pass log messages to the static member \em out and associate an intendation level and a priority with them. The implementation of MessageInterface is responsible for displaying them accordingly.
      
      Here is an example on how MessageInterface can be used. Note that for this function to work properly an implementation of MessageInterface has to be included in the final application:
  \code
  #include <core/MessageInterface.hpp> 
  
  using imaging;
  
  void some_function()
  {
    MessageInterface::out("Enter some function.", MessageInterface::IMPORTANT);
    
    MessageInterface::out("Process data.", MessageInterface::LESS_IMPORTANT, +1);
    MessageInterface::out("Data seems okay.", MessageInterface::DEBUG_ONLY, +1);
    
    MessageInterface::out("Exit some function.", MessageInterface::IMPORTANT);
  }
  \endcode
  
  */
  class MessageInterface
  {

   public:
    /** Global log message output object. */
    static MessageInterface & out;
   
    /** Different priority levels which determine how important a message is. */
    enum priority_levels {
    DEBUG_ONLY /** Debug message. */, 
    LESS_IMPORTANT /** Not so important message. */, 
    IMPORTANT /** Important message. */
    };

    virtual ~MessageInterface() {}
    
    /** Prints \em message. The priority and the level of intendation relative to the global intendation are determined by \em priority_level and \em intend. The specified intendentation is valid for this message only does not change the global intendation. */
    virtual void operator()(const std::string & message, const int priority_level, int intend = 0) = 0;

    /** Changes the global intendentation by \em intend. */
    virtual void operator()(int intend) = 0;
  };

}

#endif
