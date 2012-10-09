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

#ifndef CORE_EXCEPTION_H
#define CORE_EXCEPTION_H

#include <iostream>

namespace imaging
{
  /** \ingroup core
      \brief Generic exception with custom error message. 
      
      Base class of all exceptions thrown in \em imaging2. The user can provide a custom error message upon construction of an Exception object.
  */
  class Exception
  {
    std::string _error_msg;

  public:
    /** Constructor. The user can pass a custom error message. */
    Exception(std::string msg = std::string("Exception: An exception has ocurred!")) : _error_msg(msg) {}

    /** Returns the error message. */
    const std::string & error_msg() const { return _error_msg; }
  };

  /** \ingroup core
      \brief Mathematical exception. 
      
      An exception which is thrown because of mathematical errors, e.g. attempted division by zero.
  */
  class MathException : public Exception
  {
  public:
    /** Constructor. The user can pass a custom error message. */
    MathException(std::string msg = std::string("MathException: An exception has ocurred!")) :
    Exception(msg) {}
  };

  /** \ingroup core
      \brief File I/O exception. 
      
      An exception which is thrown because of a file related error. Most prominently these are failed attempts to open a file for reading or writing.
  */
  class FileIoException : public Exception
  {
  public:
    /** Constructor. The user can pass a custom error message. */
    FileIoException(std::string msg = std::string("FileIoException: An exception has ocurred!")) :
    Exception(msg) {}
  };
}

#endif
