// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {

  // Calculate first order derivative, dy/dt, as it appears in the differential equation
  double AlphaS_ODE::_derivative(double t, double y, const std::vector<double>& beta) const {
    /// @todo Respect the _order member variable
    const double b0 = beta[0];
    const double b1 = beta[1];
    const double b2 = beta[2];
    const double b3 = beta[3];
    const double d = (b0*y*y) + (b1*y*y*y) + (b2*y*y*y*y) + (b3*y*y*y*y*y);
    return - d / t;
  }

  // Calculate the next step, using recursion to achieve
  // adaptive step size
  vector<double> AlphaS_ODE::_rk4(double t, double y, double h, const double allowed_change) const {
    // Get beta coefficients
    const vector<double> bs = _betas(nf_Q2(t));

    // Determine increments in y based on the slopes of the function at the
    // beginning, midpoint, and end of the interval
    const double k1 = h * _derivative(t, y, bs);
    const double k2 = h * _derivative(t + h/2.0, y + k1/2.0, bs);
    const double k3 = h * _derivative(t + h/2.0, y + k2/2.0, bs);
    const double k4 = h * _derivative(t + h, y + k3, bs);
    const double next = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

    // Only call recursively if Q2 > 1 GeV (constant step
    // size after this)
    if(fabs(next) > allowed_change && t > 1.) return _rk4(t, y, h/2., allowed_change);

    vector<double> ret;
    ret.push_back(y + next);
    ret.push_back(t + h);
    return ret;
  }


  // Solve the differential equation in alphaS using an implementation of RK4
  double AlphaS_ODE::alphasQ2(double q2) const {

    /// @todo Make these class members
    // Initial step size
    /// @todo Use an adaptive step size -- stepping in log space? It needs to be *much* faster!
    /// @todo Stepping also needs to get much smaller as we approach LambdaQCD
    double h = 1;
    const double allowed_relative = 0.05;
    // Fractional threshold to which we run in Q2 scale
    /// @todo Need a mechanism to shrink the steps if accuracy < stepsize?
    /// @todo Is the accuracy in Q or Q2? Q is better, I think.
    const double faccuracy = 0.01;

    // Calculated absolute tolerance around the target Q2.
    /// @todo Would specifying the accuracy in abs rather than rel terms be better?
    double accuracy = faccuracy;

    /// @todo Running (i.e. t and h) is in Q2 or Q? Check

    /// @todo Use caching of the last target so we don't always have to start the running from MZ?
    ///       Then choose to start from last or MZ, whichever is closer.

    // Run in Q2 using RK4 algorithm until we are within our defined accuracy
    double t = sqr(_mz); // starting point
    double y = _alphas_mz; // starting value

    while (fabs(q2 - t) > accuracy) {

      /// Make the allowed change smaller as the q2 scale gets greater
      const double allowed_change = allowed_relative / sqrt(t);

      /// Mechanism to shrink the steps if accuracy < stepsize and we are close to Q2
      if (fabs(h) > accuracy && fabs(q2 - t)/h < 10 && t > 1.) h = accuracy/2.1;
      /// Take constant steps for Q2 < 1 GeV
      if (fabs(h) > 0.01 && t < 1.) {accuracy = 0.0051; h = 0.01;}

      // Check if we are going to run forward or backwards in energy scale towards target Q2.
      /// @todo C++11's std::copysign would be perfect here
      if ((q2 < t && sgn(h) > 0) || (q2 > t && sgn(h) < 0)) h *= -1;

      // Calculate next step
      vector<double> step = _rk4(t, y, h, allowed_change);

      y = step[0];
      if(y > 2.) return std::numeric_limits<double>::max();
      t = step[1];
    }
    return y;
  }

}
