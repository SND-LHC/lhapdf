// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Calculate alpha_s(Q2) by an analytic approximation
  double AlphaS_Analytic::alphasQ2(double q2) const {
    /// Get the number of active flavours and corresponding LambdaQCD
    /// @todo support any number of active flavours?
    /// Only supports active nf=3-5 at the moment

    const int nf = nf_Q2(q2);
    const double lambdaQCD = _lambdaQCD(nf);

    if(q2 <= lambdaQCD * lambdaQCD){return std::numeric_limits<double>::max();}

    // Get beta coeffs for the number of active (above threshold) quark flavours at energy Q
    const std::vector<double> beta = _betas(nf);
    const double beta02 = sqr(beta[0]);
    const double beta12 = sqr(beta[1]);

    // Pre-calculate ln(Q2/lambdaQCD) and expansion term y = 1/ln(Q2/lambdaQCD)
    const double x = q2 / (lambdaQCD*lambdaQCD);
    const double lnx = log(x);
    const double lnlnx = log(lnx);
    const double lnlnx2 = lnlnx * lnlnx;
    const double y = 1 / lnx;

    // Calculate and combine expansion terms Ay[a_0*y^0 - a_1*y^1 + By^2[a_20 + a_21 - a_22]]
    // where a_ij are series coefficients. Avoid calculating more than is needed, although it's a bit messy.
    const double A = 1 / beta[0];
    const double a_0 = 1;
    double tmp = a_0;
    if (qcdorder > 0) {
      const double a_1 = beta[1] * lnlnx / beta02;
      tmp -= a_1 * y;
    }
    if (qcdorder > 1) {
      const double B = beta12 / (beta02 * beta02);
      const double a_20 = lnlnx2 - lnlnx;
      const double a_21 = beta[2] * beta[0] / beta12;
      const double a_22 = 1;
      tmp += B * y*y * (a_20 + a_21 - a_22);
    }
    const double alphaS = A * y * tmp;
    return alphaS;
  }


}
