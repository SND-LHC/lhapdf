// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  void AlphaS_Ipol::setQ2Values(const std::vector<double>& q2s) {
    _q2s = q2s;
  }


  void AlphaS_Ipol::setAlphaSValues(const std::vector<double>& as) {
    _as = as;
  }


  /// @note This is const so it can be called silently from a const method
  void AlphaS_Ipol::_setup_grids() const {
    if (!_knotarrays.empty())
      throw LogicError("AlphaS interpolation subgrids being initialized a second time!");

    if (_q2s.size() != _as.size())
      throw MetadataError("AlphaS value and Q interpolation arrays are differently sized");

    // Walk along the Q2 vector, making subgrids at each boundary
    double blockstartQ2 = _q2s.front();
    double prevQ2 = _q2s.front();
    vector<double> q2s, as;
    for (size_t i = 0; i < _q2s.size(); ++i) {
      const double currQ2 = _q2s[i];
      const double currAS = _as[i];
      // If the Q2 value is repeated, sync the current subgrid and start a new one
      if (currQ2 == prevQ2 || i == _q2s.size()-1) {
        if (i != 0) _knotarrays[blockstartQ2] = AlphaSArray(q2s, as);
        blockstartQ2 = currQ2;
        q2s.clear(); q2s.reserve(_q2s.size() - i);
        as.clear(); as.reserve(_q2s.size() - i);
      }
      q2s.push_back(currQ2);
      as.push_back(currAS);
      prevQ2 = currQ2;
    }
  }


  double AlphaS_Ipol::_interpolateCubic(double T, double VL, double VDL, double VH, double VDH) const {
    // Pre-calculate powers of T
    const double t2 = T*T;
    const double t3 = t2*T;

    // Calculate left point
    const double p0 = (2*t3 - 3*t2 + 1)*VL;
    const double m0 = (t3 - 2*t2 + T)*VDL;

    // Calculate right point
    const double p1 = (-2*t3 + 3*t2)*VH;
    const double m1 = (t3 - t2)*VDH;

    /// @todo For consistency with the ODE solver define divergence at alpha_s > 2.
    /// @todo Will this break anything (are there PDFs that will want to interpolate values over 2. ?)
    return abs(p0 + m0 + p1 + m1) < 2. ? p0 + m0 + p1 + m1 : std::numeric_limits<double>::max();
  }


  // Interpolate alpha_s from tabulated points in Q2 via metadata
  double AlphaS_Ipol::alphasQ2(double q2) const {
    assert(q2 >= 0);

    // Use a basic constant extrapolation in case we go out of range
    if (q2 <= _q2s.front()) return _as.front(); //< @todo Gradient xpol w.r.t. logQ2 would be better, perhaps
    if (q2 >= _q2s.back()) return _as.back();

    // If this is the first valid query, set up the ipol grids
    if (_knotarrays.empty()) _setup_grids();

    // Retrieve the appropriate subgrid
    map<double, AlphaSArray>::const_iterator it = --(_knotarrays.upper_bound(q2));
    const AlphaSArray& arr = it->second;

    // Get the Q/alpha_s index on this array which is *below* this Q point
    const size_t i = arr.iq2below(q2);

    // Calculate derivatives
    double didlogq2, di1dlogq2;
    if ( i == 0 ) {
      didlogq2 = arr.ddlogq_forward(i);
      di1dlogq2 = arr.ddlogq_central(i+1);
    } else if ( i == arr.logq2s().size()-2 ) {
      didlogq2 = arr.ddlogq_central(i);
      di1dlogq2 = arr.ddlogq_backward(i+1);
    } else {
      didlogq2 = arr.ddlogq_central(i);
      di1dlogq2 = arr.ddlogq_central(i+1);
    }

    // Calculate alpha_s
    const double dlogq2 = arr.logq2(i+1) - arr.logq2(i);
    const double tlogq2 = (log(q2) - arr.logq2(i)) / dlogq2;
    return _interpolateCubic( tlogq2,
                              arr.alpha(i), didlogq2*dlogq2,
                              arr.alpha(i+1), di1dlogq2*dlogq2 );
  }


}
