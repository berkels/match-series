#ifndef __QMEXCEPTION_H
#define __QMEXCEPTION_H

// cannot include aol.h here

#include <string>
#include <iostream>
#include <sstream>

namespace {
// simple conversion utility
inline std::string tostring ( int i ) {
   std::ostringstream str;
   str << i;
   return str.str ();
 }
}

namespace aol {

/** Class for error handling. We don't want error messages to be printed if the
 *  exception is caught and consumed. Therefore we perform the dumping in the
 *  destructor and not in the constructor. Thus, we make sure that the error is
 *  printed at least once. Since we mark copies, from all copies only one
 *  (namely the last existing copy) will be dumped in its destructor.
 *
 *  This does not work correctly when throwing exceptions from within a catch block
 *  (which is perfectly ok since the old exception is considered handled upon
 *  entering the catch block). The static variables will become messed up,
 *  which may e.g. lead to printing a consumed exception in its destructor,
 *  if one throws a new one after consuming the old one.
 *
 *  \author Preusser, Droske
 */
class Exception {
public:
  /** Exception with error message Message
   */
  Exception ( const std::string& Message )
    : message ( Message ), where ( "<unknown>" ), type ( "error" ) {
    defaultInitialize();
  }

  /** Exception with error message Message and where it occured
   */
  Exception ( const std::string& Message, const std::string& Where )
    : message ( Message ), where ( Where ), type ( "error" ) {
    defaultInitialize();
  }

  /** Exception with error message Message and file and line number where
   *  it occured
   */
  Exception ( const std::string& Message, const std::string& file, int line )
    : message ( Message ), where ( "file " + file + ", line " + tostring ( line ) ), type ( "error" ) {
    defaultInitialize();
  }

  /** Exception with error message Message and file, line number and method name where
   *  it occured
   */
    Exception ( const std::string& Message, const std::string& file, int line, const std::string& methodName )
    : message ( Message ), where ( "file " + file + ", line " + tostring ( line ) + ", in method " + methodName ), type ( "error" ) {
    defaultInitialize();
  }

  /** Copy constructor marks that copies exist now.
   */
  Exception ( const Exception &other )
    : message ( other.message ), where ( other.where ), type ( other.type ) {
    wasCopied = true;
    ++copyNum;
  }

  /** Call this method to consume the exception. That is, the error will not be
   *  dumped to stderr, if the message is destructed
   */
  void consume() {
    consumed = true;
  }

  /** Destructor dumps the error if it has not been consumed yet
   */
  ~Exception() {
    if ( !wasCopied || ( wasCopied && !copyNum ) )
      if ( !consumed ) dump();
    --copyNum;
    if (!copyNum) wasCopied = false;
  }

  /** Consumes and prints the error to stderr
   */
  void dump() {
    consume();
    std::cerr << message << " : " << where << std::endl << std::flush;
  }

  /** return message string
   */
  const std::string& getMessage() {
    return message;
  }

  /** return exception type
   */
  const std::string getType() {
    return type;
  }

  /** return, where the exception occured
   */
  const std::string getWhere() {
    return where;
  }

protected:
  /** Derived classes use this constructor to further specify their type
   */
  Exception ( const std::string& myName, const std::string& method, const std::string& file, int line )
    : message ( myName + " in " + method ), where ( "file " + file + ", line " + tostring ( line ) ), type ( "error" ) {
    consumed = wasCopied = false;
    copyNum = 0;
  }

private:
  void defaultInitialize ( ) {
    consumed = wasCopied = false;
    copyNum = 0;
  }

protected:
  //! The message of this exception
  std::string message;
  //! Where the exception occured
  std::string where;
  //! Type of exception
  std::string type;

private:
  static bool consumed, wasCopied;    // Defined in Quoc.cpp
  static int  copyNum;                // Defined in Quoc.cpp
};

// Input and output

class IOException : public Exception {
public:
  IOException ( const std::string& method, const std::string& file, int line ) : Exception ( "IOException", method, file, line ) {}};

class PipeException : public Exception {
public:
  PipeException ( const std::string& method, const std::string& file, int line ) : Exception ( "PipeException", method, file, line ) {}};

class FileException : public Exception {
public:
  FileException ( const std::string& method, const std::string& file, int line ) : Exception ( "FileException", method, file, line ) {}};

class FileFormatException : public Exception {
public:
  FileFormatException ( const std::string& method, const std::string& file, int line ) : Exception ( "FileFormatException", method, file, line ) {}};

class FileDoesNotExistException : public Exception {
public:
  FileDoesNotExistException ( const std::string& method, const std::string& file, int line ) : Exception ( "FileDoesNotExistException", method, file, line ) {}};

class DicomImportException : public Exception {
public:
  DicomImportException ( const std::string& method, const std::string& file, int line ) : Exception ( "DicomImportException", method, file, line ) {}};

// Bad parameters

class TypeException : public Exception {
public:
  TypeException ( const std::string& method, const std::string& file, int line ) : Exception ( "TypeException", method, file, line ) {}};

class InconsistentDataException : public Exception {
public:
  InconsistentDataException ( const std::string& method, const std::string& file, int line ) : Exception ( "InconsistentDataException", method, file, line ) {}};

class ParameterException : public Exception {
public:
  ParameterException ( const std::string& method, const std::string& file, int line ) : Exception ( "ParameterException", method, file, line ) {}};

class DimensionMismatchException : public Exception {
public:
  DimensionMismatchException ( const std::string& method, const std::string& file, int line ) : Exception ( "DimensionMismatchException", method, file, line ) {}};

class ParameterRangeException : public Exception {
public:
  ParameterRangeException ( const std::string& method, const std::string& file, int line ) : Exception ( "ParameterRangeException", method, file, line ) {}};

class OutOfBoundsException : public Exception {
public:
  OutOfBoundsException ( const std::string& method, const std::string& file, int line ) : Exception ( "OutOfBoundsException", method, file, line ) {}};

// Runtime errors

class OutOfMemoryException : public Exception {
public:
  OutOfMemoryException ( const std::string& method, const std::string& file, int line ) : Exception ( "OutOfMemoryException", method, file, line ) {}};

class UnimplementedCodeException : public Exception {
public:
  UnimplementedCodeException ( const std::string& method, const std::string& file, int line ) : Exception ( "UnimplementedCodeException", method, file, line ) {}};

// Logical errors

class OperationNotPermittedException : public Exception {
public:
  OperationNotPermittedException ( const std::string& method, const std::string& file, int line ) : Exception ( "OperationNotPermittedException", method, file, line ) {}};

class LogicException : public Exception {
public:
  LogicException ( const std::string& method, const std::string& file, int line ) : Exception ( "LogicException", method, file, line ) {}};


/**
 * Custom terminate handler intended to catch unhandled exceptions. In particular useful
 * to get hold of exceptions thrown inside OpenMP blocks.
 *
 *  \author Berkels
 */
void catchGlobal();

} // end namespace

#endif // __QMEXCEPTION_H
