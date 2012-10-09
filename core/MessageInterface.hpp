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
