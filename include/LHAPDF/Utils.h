#pragma once

// STL includes
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
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


  // String handling
  template <typename T>
  inline std::string to_str(const T& val) {
    return lexical_cast<string>(val);
  }
  template <typename T>
  inline std::string to_str(const std::vector<T>& vec) {
    vector<string> svec; svec.reserve(vec.size());
    foreach (const T& t, vec) svec.push_back( to_str(t) );
    return join(svec, ",");
  }
  inline std::string to_str_zeropad(int val, size_t nchars=4) {
    stringstream ss;
    ss << setfill('0') << setw(nchars) << val;
    return ss.str();
  }


  // Math functions in the LHAPDF namespace

  /// Convenience function for squaring (of any type)
  template <typename N>
  inline N sqr(const N& x) { return x*x; }

  /// @todo Add iszero() & equals(,) functions?


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


}
