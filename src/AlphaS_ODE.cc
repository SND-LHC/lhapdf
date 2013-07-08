// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  /// @todo Respect the NumFlavors + FlavorScheme metadata


  // Calculate first order derivative, dy/dt, as it appears in the differential equation
  double AlphaS_ODE::_derivative(double t, double y, const vector<double>& beta) const {
    /// @todo Respect the _order member variable
    const double b0 = beta[0];
    const double b1 = beta[1];
    const double b2 = beta[2];
    const double b3 = beta[3];
    const double d = (b0*y*y) + (b1*y*y*y) + (b2*y*y*y*y) + (b3*y*y*y*y*y);
    return - d / t;
  }


  // Calculate the next step, using recursion to achieve
  // adaptive step size. Passing be reference explained
  // below.
  void AlphaS_ODE::_rk4(double& t, double& y, double h,
                        const double allowed_change, const vector<double>& bs) const {

    // Determine increments in y based on the slopes of the function at the
    // beginning, midpoint, and end of the interval
    const double k1 = h * _derivative(t, y, bs);
    const double k2 = h * _derivative(t + h/2.0, y + k1/2.0, bs);
    const double k3 = h * _derivative(t + h/2.0, y + k2/2.0, bs);
    const double k4 = h * _derivative(t + h, y + k3, bs);
    const double change = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

    // Only call recursively if Q2 > 1 GeV (constant step
    // size after this)
    if (t > 1. && fabs(change) > allowed_change) {
      _rk4(t, y, h/2., allowed_change, bs);
    } else {
      y += change;
      t += h;
    }
  }


  // Solve for alpha_s(q2) using alpha_s = y and Q2 = t as starting points
  // Pass y and t by reference and change them to avoid having to
  // return them in a container -- bit confusing but should be more
  // efficient
  void AlphaS_ODE::_solve(double q2, double& t, double& y,
                          const double& allowed_relative, double h, double accuracy) const {
    while (fabs(q2 - t) > accuracy) {
      /// Make the allowed change smaller as the q2 scale gets greater
      const double allowed_change = allowed_relative / t;

      /// Mechanism to shrink the steps if accuracy < stepsize and we are close to Q2
      if (fabs(h) > accuracy && fabs(q2 - t)/h < 10 && t > 1.) h = accuracy/2.1;
      /// Take constant steps for Q2 < 1 GeV
      if (fabs(h) > 0.01 && t < 1.) {accuracy = 0.00051; h = 0.001;}
      // Check if we are going to run forward or backwards in energy scale towards target Q2.
      /// @todo C++11's std::copysign would be perfect here
      if ((q2 < t && sgn(h) > 0) || (q2 > t && sgn(h) < 0)) h *= -1;

      // Get beta coefficients
      const vector<double> bs = _betas(numFlavorsQ2(t));

      // Calculate next step
      _rk4(t, y, h, allowed_change, bs);

      if (y > 20.) { y = std::numeric_limits<double>::max(); break; }
    }
  }


  /// Interpolate to get Alpha_S if the ODE has been solved,
  /// otherwise solve ODE from scratch
  double AlphaS_ODE::alphasQ2(double q2) const {
    // Tabulate ODE solutions for interpolation and return interpolated value
    _interpolate();
    return _ipol.alphasQ2(q2);
    // // Or directly return the ODE result (for testing)
    // double h = 2;
    // const double allowed_relative = 0.01;
    // const double faccuracy = 0.01;
    // double accuracy = faccuracy;
    // double t = sqr(_mz); // starting point
    // double y = _alphas_mz; // starting value
    // _solve(q2, t, y, allowed_relative, h, accuracy); // solve ODE
    // return y;
  }


  // Solve the differential equation in alphaS using an implementation of RK4
  //
  /// @todo Stepping also needs to get much smaller as we approach LambdaQCD
  /// -- Can't use adaptive step size close to lambdaQCD effectively (since changes will be large)
  /// -- Do they really need to be *very* small (since errors always will get large there)?
  void AlphaS_ODE::_interpolate() const {
    if ( _calculated ) return;

    // Initial step size
    double h = 2.0;
    /// This the the relative error allowed for the adaptive step size. Should be optimised.
    const double allowed_relative = 0.01;
    // Fractional threshold to which we run in Q2 scale
    /// @todo Is the accuracy in Q or Q2?
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

    // If a vector of anchor points in q2 has been given, solve for those.
    if ( !_q2s.empty() ) {

      vector<double> alphas;
      foreach (double q2, _q2s) {
        _solve(q2, t, y, allowed_relative, h, accuracy);
        alphas.push_back(y);
        // If alpha_s goes over 20 the ODE diverges too fast, so go back to start
        // rather than using this point as a starting point for the next one
        if (y > 20.) { t = sqr(_mz); y = _alphas_mz; }
      }
      _ipol.setQ2Values(_q2s);
      _ipol.setAlphaSValues(alphas);

    } else {

      // Start evolution in Q at MZ, and assemble a grid of anchor points.
      vector<double> q2s, alphas;

      // To save time we solve from MZ down to Q=0.5, then go back to MZ and solve up to Q=1000
      double tmp = sqr(_mz); //< @todo Use a better variable name to indicate that this is the scale evolution variable
      while (tmp > sqr(0.5)) {
        _solve(tmp, t, y, allowed_relative, h, accuracy);
        q2s.push_back(t);
        alphas.push_back(y);
        if (y > 20.) break;
        tmp -= (10 * accuracy * t);
      }
      t = sqr(_mz); // starting point
      y = _alphas_mz; // starting value
      tmp = sqr(_mz);
      while (tmp < sqr(1000)) {
        tmp += (10 * accuracy * t);
        _solve(tmp, t, y, allowed_relative, h, accuracy);
        q2s.push_back(t);
        alphas.push_back(y);
      }

      // We assume alpha_s is monotonic
      std::sort(q2s.begin(), q2s.end());
      std::sort(alphas.begin(), alphas.end(), cmpDescend<double>);

      // Set interpolation knots and values
      _ipol.setQ2Values(q2s);
      _ipol.setAlphaSValues(alphas);
    }

    _calculated = true;
  }


}
