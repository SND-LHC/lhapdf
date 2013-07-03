// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {

  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS_Ipol::nf_Q2(double q2) const {
    int nf = 0;
    for (int it = 0; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return nf;
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

    return p0 + m0 + p1 + m1;
  }

  double AlphaS_Ipol::_ddq_forward( size_t i ) const {
    return (_as[i+1] - _as[i]) / (_q2s[i+1] - _q2s[i]);
  }

  double AlphaS_Ipol::_ddq_backward( size_t i ) const {
    return (_as[i] - _as[i-1]) / (_q2s[i] - _q2s[i-1]);
  }

  double AlphaS_Ipol::_ddq_central( size_t i ) const {
    return (_ddq_forward(i) + _ddq_backward(i))*0.5;
  }

  // Interpolate alpha_s from tabulated points in Q2 via metadata
  double AlphaS_Ipol::alphasQ2(double q2) const {
    //Cubic interpolation of std::vector<double> q2 and as
    size_t i = std::upper_bound( _q2s.begin(), _q2s.end(), q2 ) - _q2s.begin();
    if( i == _q2s.size()) i -= 1;
    i -= 1;

    //Calculate derivatives
    double didq2, di1dq2;
    if( i == 0 ) {
        didq2 = _ddq_forward(i);
        di1dq2 = _ddq_central(i+1);
    }
    else if( i == _q2s.size()-1 ) {
        didq2 = _ddq_central(i);
        di1dq2 = _ddq_backward(i+1);
    }
    else {
        didq2 = _ddq_central(i);
        di1dq2 = _ddq_central(i+1);
    }
    
    //Calculate AlphaS
    double dq2 = _q2s[i+1] - _q2s[i];
    double tq2 = (q2 - _q2s[i]) / dq2;
    return _interpolateCubic( tq2, _q2s[i], didq2, _q2s[i+1], di1dq2 );
  }
}
