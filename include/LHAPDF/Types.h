#pragma once

#include <stdint.h> // -> cstdint with C++11
#include <cstdlib>
#include <vector>
#include <map>

#include "boost/foreach.hpp"
#ifndef foreach
#define foreach BOOST_FOREACH
#endif

namespace LHAPDF {

  /// Safe type for PID codes (is int safe enough?)
  typedef long PID_t;

  /// @todo Kill
  typedef double* PIDdata;

  /// @todo Kill
  typedef std::map<PID_t,PIDdata> PIDmap;

  /// @todo Kill
  typedef std::vector<double> AxisKnots;

}
