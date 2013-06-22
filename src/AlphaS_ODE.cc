// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {

  void AlphaS_ODE::setLambda(unsigned int i, double lambda){throw Exception ("Can't set lambdas for ODE solver.");}

  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS_ODE::nf_Q2(double q2) const {
    int nf = 0;
    for (int it = 0; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return nf;
  }

  namespace { // unnamed namespace

    // Calculate first order derivative, dy/dt, as it appears in the differential equation
    double _derivative(double t, double y, const std::vector<double>& beta) {
      /// @todo Respect the _order member variable
      const double d0 = beta[0]/(-2*M_PI);
      const double d1 = beta[1]/(-4*M_PI*M_PI);
      const double d2 = beta[2]/(-64*M_PI*M_PI*M_PI);
      const double d = (d0*y*y) + (d1*y*y*y) + (d2*y*y*y*y);
      return d / t;
    }

  }


  // Solve the differential equation in alphaS using an implementation of RK4
  double AlphaS_ODE::alphasQ2(double q2) const {

    /// @todo Make these class members
    // Initial step size
    /// @todo Use an adaptive step size -- stepping in log space? It needs to be *much* faster!
    /// @todo Stepping also needs to get much smaller as we approach LambdaQCD
    double h = 0.01;
    // Fractional threshold to which we run in Q2 scale
    /// @todo Need a mechanism to shrink the steps if accuracy < stepsize?
    /// @todo Is the accuracy in Q or Q2? Q is better, I think.
    const double faccuracy = 0.001;

    // Calculated absolute tolerance around the target Q2.
    /// @todo Would specifying the accuracy in abs rather than rel terms be better?
    const double accuracy = faccuracy * q2;

    /// @todo Running (i.e. t and h) is in Q2 or Q? Check

    /// @todo Use caching of the last target so we don't always have to start the running from MZ?
    ///       Then choose to start from last or MZ, whichever is closer.

    // Run in Q2 using RK4 algorithm until we are within our defined accuracy
    double t = sqr(_mz); // starting point
    double y = _alphas_mz; // starting value
    while (fabs(q2 - t) > accuracy) {
      // Calculate beta functions based on the number of active flavours at Q2=t
      const vector<double> bs = _betas(nf_Q2(t));

      /// Mechanism to shrink the steps if accuracy < stepsize and we are close to Q2
      if (fabs(h) > accuracy && fabs(q2 - t)/h < 10) h = accuracy/2.0;

      // Check if we are going to run forward or backwards in energy scale towards target Q2.
      /// @todo C++11's std::copysign would be perfect here
      if ((q2 < t && sgn(h) > 0) || (q2 > t && sgn(h) < 0)) h *= -1;

      // cout << t << " -> " << q2 << ", in steps of size " << h << endl;

      // Determine increments in y based on the slopes of the function at the
      // beginning, midpoint, and end of the interval
      const double k1 = h * _derivative(t, y, bs);
      const double k2 = h * _derivative(t + h/2.0, y + k1/2.0, bs);
      const double k3 = h * _derivative(t + h/2.0, y + k2/2.0, bs);
      const double k4 = h * _derivative(t + h, y + k3, bs);

      // Calculate the weighted average of these increments and make the step
      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      t += h;
    }
    return y;
  }


}
