// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


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
