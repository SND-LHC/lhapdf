#pragma once

#include <exception>
#include <stdexcept>

namespace LHAPDF {


  /// @name Exception classes for error handling
  //@{

  /// @brief Generic unspecialised LHAPDF runtime error.
  ///
  /// NB. We don't use "Error" because that has a physics meaning!
  class Exception : public std::runtime_error {
  public:
    Exception(const std::string& what) : std::runtime_error(what) {}
  };


  /// Error for general PDF grid problems.
  class GridError : public Exception {
  public:
    GridError(const std::string& what) : Exception(what) {}
  };


  /// Error to be thrown when out of the valid range of a PDF.
  class RangeError : public Exception {
  public:
    RangeError(const std::string& what) : Exception(what) {}
  };


  /// Error for places where it should not have been possible to get to!
  class LogicError : public Exception {
  public:
    LogicError(const std::string& what) : Exception(what) {}
  };


  /// @brief Error for unfound or broken metadata entries.
  class MetadataError : public Exception {
  public:
    MetadataError(const std::string& what) : Exception(what) {}
  };


  /// @brief Error for file reading errors.
  class ReadError : public Exception {
  public:
    ReadError(const std::string& what) : Exception(what) {}
  };


  /// @brief Error for requests for unsupported/invalid flavour PIDs.
  class FlavorError : public Exception {
  public:
    FlavorError(const std::string& what) : Exception(what) {}
  };


  /// @brief Error to be raised by object factories given invalid requests.
  class FactoryError : public Exception {
  public:
    FactoryError(const std::string& what) : Exception(what) {}
  };


  /// @brief Error to be raised when a LHAPDF ID indexing fails
  class IndexError : public Exception {
  public:
    IndexError(const std::string& what) : Exception(what) {}
  };


  /// Problem exists between keyboard and chair.
  class UserError : public Exception {
  public:
    UserError(const std::string& what) : Exception(what) {}
  };

  //@}


}
