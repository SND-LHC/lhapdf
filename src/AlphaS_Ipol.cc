// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  /// @todo Add subgrid support


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

    return p0 + m0 + p1 + m1;
  }


  double AlphaS_Ipol::_ddq_forward( size_t i ) const {
    return (_as[i+1] - _as[i]) / (_logq2s[i+1] - _logq2s[i]);
  }

  double AlphaS_Ipol::_ddq_backward( size_t i ) const {
    return (_as[i] - _as[i-1]) / (_logq2s[i] - _logq2s[i-1]);
  }

  double AlphaS_Ipol::_ddq_central( size_t i ) const {
    return (_ddq_forward(i) + _ddq_backward(i))*0.5;
  }


  // Interpolate alpha_s from tabulated points in Q2 via metadata
  double AlphaS_Ipol::alphasQ2(double q2) const {
    // Actually we operate in log space
    const double logq2 = log(q2);

    // Make sure the values to interpolate from have the same dimensions
    assert(_as.size() == _logq2s.size());

    // Use a basic constant extrapolation in case we go out of range
    if (logq2 <= _logq2s.front()) return _as.front();
    if (logq2 >= _logq2s.back()) return _as.back();

    // Cubic interpolation of std::vector<double> q2 and as
    size_t i = std::upper_bound( _logq2s.begin(), _logq2s.end(), logq2 ) - _logq2s.begin();
    if (i == _logq2s.size()) i -= 1;
    i -= 1;

    // Calculate derivatives
    double didlogq2, di1dlogq2;
    if ( i == 0 ) {
      didlogq2 = _ddq_forward(i);
      di1dlogq2 = _ddq_central(i+1);
    } else if ( i == _logq2s.size()-2 ) {
      didlogq2 = _ddq_central(i);
      di1dlogq2 = _ddq_backward(i+1);
    } else {
      didlogq2 = _ddq_central(i);
      di1dlogq2 = _ddq_central(i+1);
    }

    // Calculate alpha_s
    double dlogq2 = _logq2s[i+1] - _logq2s[i];
    double tlogq2 = (logq2 - _logq2s[i]) / dlogq2;
    return _interpolateCubic( tlogq2, _as[i], didlogq2, _as[i+1], di1dlogq2 );
  }


}
