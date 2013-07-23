// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

// STL includes
#include <cassert>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
// Boost includes
#include "boost/shared_ptr.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/range.hpp"
// #include "boost/range/algorithm/sort.hpp"
// #include "boost/range/algorithm/lower_bound.hpp"
// #include "boost/range/algorithm/upper_bound.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"
#include "boost/bind.hpp"

// Cosmetic wrapper for BOOST_FOREACH macro (until we can use C++11 range-based for loops)
#ifndef foreach
#define foreach BOOST_FOREACH
#endif


/// Namespace for all LHAPDF functions and classes
namespace LHAPDF {
  // Allow implicit use of the std and boost namespaces within namespace LHAPDF
  using namespace std;
  using namespace boost;


  /// @name String handling utility functions
  //@{

  /// Make a string representation of @a val
  template <typename T>
  inline std::string to_str(const T& val) {
    return lexical_cast<string>(val);
  }

  /// Make a string representation of a vector @a vec
  template <typename T>
  inline std::string to_str(const std::vector<T>& vec) {
    vector<string> svec; svec.reserve(vec.size());
    foreach (const T& t, vec) svec.push_back( to_str(t) );
    return join(svec, ",");
  }

  /// Format an integer @a val as a zero-padded string of length @a nchars
  inline std::string to_str_zeropad(int val, size_t nchars=4) {
    stringstream ss;
    ss << setfill('0') << setw(nchars) << val;
    return ss.str();
  }

  //@}


  /// @name Math functions in the LHAPDF namespace
  //@{

  /// Convenience function for squaring (of any type)
  template <typename N>
  inline N sqr(const N& x) { return x*x; }

  /// Get the sign of a number
  template <typename N>
  inline int sgn(N val) { return (N(0) < val) - (val < N(0)); }

  /// @todo Add iszero() & equals(,) functions?

  //@}

  /// @name Container handling helpers
  //@{

  // /// @name Implementation of generic begin/end container identification by traits
  // /// taken from http://stackoverflow.com/a/9407420/91808 . Needs C++11 (or maybe just C++0x).
  // //@{

  // #include <type_traits>
  // //#include "boost/type_traits.hpp"

  // template<typename T>
  // struct has_const_iterator {
  // private:
  //   typedef char                      yes;
  //   typedef struct { char array[2]; } no;
  //   template<typename C> static yes test(typename C::const_iterator*);
  //   template<typename C> static no  test(...);
  // public:
  //   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
  //   typedef T type;
  // };

  // template <typename T>
  // struct has_begin_end {
  //   template<typename C> static char (&f(typename std::enable_if<
  //     std::is_same<decltype(static_cast<typename C::const_iterator (C::*)() const>(&C::begin)),
  //     typename C::const_iterator(C::*)() const>::value, void>::type*))[1];

  //   template<typename C> static char (&f(...))[2];

  //   template<typename C> static char (&g(typename std::enable_if<
  //     std::is_same<decltype(static_cast<typename C::const_iterator (C::*)() const>(&C::end)),
  //     typename C::const_iterator(C::*)() const>::value, void>::type*))[1];

  //   template<typename C> static char (&g(...))[2];

  //   static bool const beg_value = sizeof(f<T>(0)) == 1;
  //   static bool const end_value = sizeof(g<T>(0)) == 1;
  // };

  // template<typename T>
  // struct is_container
  //   : std::integral_constant<bool, has_const_iterator<T>::value &&
  //                            has_begin_end<T>::beg_value && has_begin_end<T>::end_value>
  // { };

  //@}

}
