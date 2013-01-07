#pragma once

// STL includes
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
// Boost includes
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"

// Cosmetic wrapper for BOOST_FOREACH macro (until we can use C++11 range-based for loops)
#ifndef foreach
#define foreach BOOST_FOREACH
#endif


namespace LHAPDF {
  // Allow implicit use of the std and boost namespaces within namespace LHAPDF
  using namespace std;
  using namespace boost;


  // Boost string handling
  template <typename T>
  std::string to_str(const T& val) {
    return lexical_cast<std::string>(val);
  }


  // Math functions in the LHAPDF namespace

  /// Convenience function for squaring (of any type)
  template <typename N>
  inline N sqr(const N& x) { return x*x; }

  /// @todo Add iszero() & equals(,) functions?
}
