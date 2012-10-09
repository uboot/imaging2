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

#ifndef CORE_CMESSAGE_H
#define CORE_CMESSAGE_H

#include <core/MessageInterface.hpp>
#include <iostream>

namespace imaging
{

  /** \ingroup core
      \brief Implementation of MessageInterface.
       
      This class provides an implementation of the \em imaging2 log message interface. The class defines the static member <em>CMessage::out</em> is declared. Cmessage formats them accordingly and prints them to the standard output.
      
      The user can adjust the verbose level. This function is specific to Cmessage and not part of MessageInterface. Ideally the user includes <tt>core/Cmessage.hpp</tt> in as little files as possible and mainly relies on <tt>core/MessageInterface.hpp</tt> for displaying log messages.
      
      Here is an example on how Cmessage can be used:
  \code
  #include <core/Cmessage.hpp> 
  
  using imaging;
  
  int main ( int argc, char **argv )
  {
    // configure Cmessage::cmessage_out to display ALL messages
    // this function is specific to Cmessage
    Cmessage::out.set_verbose_level(Cmessage::DEBUG_MESSAGES);
 
    // from now on only MessageInterface is required
    MessageInterface::out("Start program.", MessageInterface::IMPORTANT);
    
    // user enters data
    
    MessageInterface::out("User entered:", MessageInterface::LESS_IMPORTANT, +1);
    MessageInterface::out(data, MessageInterface::LESS_IMPORTANT, +1);
    
    // check if data is valid
    
    MessageInterface::out(+1)
    MessageInterface::out("Wrong data entered!", MessageInterface::DEBUG_ONLY);
    MessageInterface::out("User is an idiot!", MessageInterface::DEBUG_ONLY);
    MessageInterface::out(-1)
    
    MessageInterface::out("Finished program.", MessageInterface::IMPORTANT);
  }
  \endcode
  
      Under the condition that \em MessageInterface::out is configured to display all messages (including those of type <em>MessageInterface::DEBUG_ONLY</em>) the above code will produce the following result:
  \code
  Start program.
    User entered:
    ...the data...
      Wrong data entered!
      User is an idiot!
  Finished program.
  \endcode
  */
  class Cmessage : public MessageInterface
  {

    int _current_intendation_level;
    size_t _current_verbose_level;

  public:
    static Cmessage out;

    /** Different verbose levels which determine which messages are actually displayed. */
    enum verbose_levels { 
      DEBUG_MESSAGES /** Display all messages including debug messages.*/,
      ALL_MESSAGES /** Display all messages but no debug messages. */, 
      LITTLE_MESSAGES /** Display only important messages. */ , 
      NO_MESSAGES /** Display no messages. */
    };

    Cmessage();

    /** Sets the current verbose level of the message object. This is e.g. useful, if the user wants to temporarily increase the verbose level. */
    void set_verbose_level(size_t level);
    
    /** Returns the current verbose level. */
    const size_t verbose_level() const { return _current_verbose_level; }
    
    void operator()(const std::string & message, const int priority_level, int intend = 0);

    void operator()(int intend);
  };

}

#endif
