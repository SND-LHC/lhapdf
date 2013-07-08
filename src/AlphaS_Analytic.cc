// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Calculate the number of active quark flavours at energy scale Q2.
  // Respects min/max nf
  int AlphaS_Analytic::numFlavorsQ2(double q2) const {
    int nf = _nfmin;
    for (int it = _nfmin; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmasses[it]) && _qmasses[it] != 0) nf = it+1;
    return (nf > _nfmax) ? _nfmax : nf;
  }


  // Set lambda_i && recalculate nfmax and nfmin
  void AlphaS_Analytic::setLambda(unsigned int i, double lambda) {
    if (_lambdas.size() < i)_lambdas.resize(i);
    _lambdas[i-1] = lambda;
    _setFlavors();
  }


  // Recalculate nfmax and nfmin after a new lambda has been set
  // Relies on vector<double> initialising all elements to 0 by
  // default, I *think* this should work even though I am comparing
  // floats since 0 has a precise representation
  void AlphaS_Analytic::_setFlavors() {
    for (size_t it = 0; it < _lambdas.size(); ++it) {
      if (_lambdas[it] != 0.) { _nfmin = it+1; break; }
    }
    for (size_t it = _lambdas.size() - 1; it >= 0; --it) {
      if (_lambdas[it] != 0.) { _nfmax = it+1; break; }
    }
    if (_nfmin != _nfmax) {
      for (int it = _nfmin; it < _nfmax - 1; ++it) {
        if (_lambdas[it] == 0) throw Exception ("Need to set intermediate lambdas.");
      }
    }
  }


  // Return the correct lambda for a given number of active flavours
  double AlphaS_Analytic::_lambdaQCD(int nf) const {
    double lambda = _lambdas[nf-1];
    if (lambda == 0)
      throw Exception("Invalid nf " + to_str(nf) + " requested for lambdaQCD");
    return lambda;
  }


  // Calculate alpha_s(Q2) by an analytic approximation
  double AlphaS_Analytic::alphasQ2(double q2) const {
    /// Get the number of active flavours and corresponding LambdaQCD
    /// @todo support any number of active flavours?
    /// Only supports active nf=3-5 at the moment

    const int nf = numFlavorsQ2(q2);
    const double lambdaQCD = _lambdaQCD(nf);

    if (q2 <= lambdaQCD * lambdaQCD) return std::numeric_limits<double>::max();

    // Get beta coeffs for the number of active (above threshold) quark flavours at energy Q
    const std::vector<double> beta = _betas(nf);
    const double beta02 = sqr(beta[0]);
    const double beta12 = sqr(beta[1]);

    // Pre-calculate ln(Q2/lambdaQCD) and expansion term y = 1/ln(Q2/lambdaQCD)
    const double x = q2 / (lambdaQCD*lambdaQCD);
    const double lnx = log(x);
    const double lnlnx = log(lnx);
    const double lnlnx2 = lnlnx * lnlnx;
    const double lnlnx3 = lnlnx * lnlnx * lnlnx;
    const double y = 1 / lnx;

    // Calculate terms up to qcdorder = 3
    // A bit messy because the actual expressions are
    // quite messy...
    const double A = 1 / beta[0];
    const double a_0 = 1;
    double tmp = a_0;
    if (_qcdorder > 0) {
      const double a_1 = beta[1] * lnlnx / beta02;
      tmp -= a_1 * y;
    }
    if (_qcdorder > 1) {
      const double B = beta12 / (beta02 * beta02);
      const double a_20 = lnlnx2 - lnlnx;
      const double a_21 = beta[2] * beta[0] / beta12;
      const double a_22 = 1;
      tmp += B * y*y * (a_20 + a_21 - a_22);
    }
    if (_qcdorder > 2) {
      const double C = 1. / (beta02 * beta02 * beta02);
      const double a_30 = (beta12 * beta[1]) * (lnlnx3 - (5/2.) * lnlnx2 - 2 * lnlnx + 0.5);
      const double a_31 = 3 * beta[0] * beta[1] * beta[2] * lnlnx;
      const double a_32 = 0.5 * beta02 * beta[3];
      tmp -= C * y*y*y * (a_30 + a_31 - a_32);
    }
    const double alphaS = A * y * tmp;
    return alphaS;
  }


}
