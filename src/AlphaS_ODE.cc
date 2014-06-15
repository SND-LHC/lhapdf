// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"
#include "boost/bind.hpp"

namespace LHAPDF {

  // Calculate first order derivative, dy/dt, as it appears in the differential equation
  double AlphaS_ODE::_derivative(double t, double y, const vector<double>& beta) const {
    if ( _qcdorder == 0 ) return 0;
    const double b0 = beta[0];
    double d = (b0*y*y);
    if ( _qcdorder == 1 ) return - d / t;
    const double b1 = beta[1];
    d += (b1*y*y*y);
    if ( _qcdorder == 2 ) return - d / t;
    const double b2 = beta[2];
    d += (b2*y*y*y*y);
    if ( _qcdorder == 3 ) return - d / t;
    const double b3 = beta[3];
    d += (b3*y*y*y*y*y);
    return - d / t;
  }


  // Calculate the next step, using recursion to achieve
  // adaptive step size. Passing by reference explained
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
    if ( q2 == t ) return;
    while (fabs(q2 - t) > accuracy) {
      /// Make the allowed change smaller as the q2 scale gets greater
      const double allowed_change = allowed_relative / t;

      /// Mechanism to shrink the steps if accuracy < stepsize and we are close to Q2
      if (fabs(h) > accuracy && fabs(q2 - t)/h < 10 && t > 1.) h = accuracy/2.1;
      /// Take constant steps for Q2 < 1 GeV
      if (fabs(h) > 0.01 && t < 1.) { accuracy = 0.00051; h = 0.001; }
      // Check if we are going to run forward or backwards in energy scale towards target Q2.
      /// @todo C++11's std::copysign would be perfect here
      if ((q2 < t && sgn(h) > 0) || (q2 > t && sgn(h) < 0)) h *= -1;

      // Get beta coefficients
      const vector<double> bs = _betas(numFlavorsQ2(t));

      // Calculate next step
      _rk4(t, y, h, allowed_change, bs);

      if (y > 2.) { y = std::numeric_limits<double>::max(); break; }
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
  void AlphaS_ODE::_interpolate() const {
    if ( _calculated ) return;

    // Initial step size
    double h = 2.0;
    /// This the the relative error allowed for the adaptive step size.
    /// @todo Should be optimised.
    const double allowed_relative = 0.01;

    /// Accuracy of Q2 (error in Q2 within this / 2)
    double accuracy = 0.01;

    // Run in Q2 using RK4 algorithm until we are within our defined accuracy
    double t = sqr(_mz); // starting point
    double y = _alphas_mz; // starting value

    // If a vector of knots in q2 has been given, solve for those.
    if ( !_q2s.empty() ) {

      vector<double> alphas;
      double low_lim = 0;
      double last_val = -1;
      //      BOOST_FOREACH (double q2, _q2s)
      for ( size_t ind = 0; ind < _q2s.size(); ++ind) {
        double q2 = _q2s[ind];
        /// Deal with cases with two identical adjacent points (thresholds) by decreasing step size,
        /// allowed errors, and accuracy.
        if ( ind != _q2s.size() - 1 ) {
          if ( q2 == _q2s[ind+1] ) {
            last_val = q2;
            _solve(q2, t, y, allowed_relative/5, h/5, accuracy/5);
            alphas.push_back(y);
            // Define divergence after y > 2. -- we have no accuracy after that any way
            if ( y > 2. ) { low_lim = q2; }
            continue;
          }
        }
        // If q2 is lower than a value that already diverged, it will also diverge
        if ( q2 < low_lim ) {
          alphas.push_back( std::numeric_limits<double>::max() );
          t = sqr(_mz);
          y = _alphas_mz;
          continue;
        // If last point was the same we don't need to recalculate
        } else if ( q2 == last_val ) {
          alphas.push_back(y);
          continue;
        // Else calculate
        } else {
          last_val = q2;
          _solve(q2, t, y, allowed_relative, h, accuracy);
          alphas.push_back(y);
          // Define divergence after y > 2. -- we have no accuracy after that any way
          if ( y > 2. ) { low_lim = q2; }
        }
      }

      _ipol.setQ2Values(_q2s);
      _ipol.setAlphaSValues(alphas);

    } else {

      // Start evolution in Q at MZ, and assemble a grid of anchor points.
      vector<double> q2s, alphas;
      vector<pair<int, double> > grid; // for storing in correct order
      int index = 0; // for sorting in the correct order

      // To save time we solve from MZ down to Q=0.5, then go back to MZ and solve up to Q=1000
      /// @todo The knots are not optimised at the moment
      double knot = sqr(_mz);
      while ( knot > sqr(0.5) ) {
        _solve(knot, t, y, allowed_relative, h, accuracy);
        q2s.push_back(t);
//        alphas.push_back(y);
        grid.push_back(make_pair(index, y));
        index--;
        if ( y > 2. ) break;
        knot -= (10 * accuracy * t);
      }
      t = sqr(_mz); // starting point
      y = _alphas_mz; // starting value
      knot = sqr(_mz);
      index = 1;
      while ( knot < sqr(1000) ) {
        knot += (10 * accuracy * t);
        _solve(knot, t, y, allowed_relative, h, accuracy);
        q2s.push_back(t);
//        alphas.push_back(y);
        grid.push_back(make_pair(index, y));
        index++;
      }

      // Sorting the values in the correct order
      /// @todo AB: Probably this "magic" isn't needed? Not that it's wrong, but it is obscure -- at least add an explanatory comment!
      /// KN: This just sorts the calculated alpha_s values into the correct order (since I interpolate first down from M_Z to sqrt(0.5), and then
      /// go back to M_Z and go up). There's probably a cleaner way to do it, this is a modified StackExchange suggestion!
      /// It sorts the vector<pair<int, double> > by ascending int values, it is required since the interpolator assumes the vectors to be ordered
      std::sort(grid.begin(), grid.end(),
                boost::bind(&std::pair<int, double>::first, _1) < boost::bind(&std::pair<int, double>::first, _2));

      // Need to do this since BOOST_FOREACH can't deal with arguments
      // with ,s (and I don't want to clutter with a typedef)
      for ( size_t x = 0; x < grid.size(); ++x ) {
         alphas.push_back(grid.at(x).second);
      }

      std::sort(q2s.begin(), q2s.end());
//      std::sort(alphas.begin(), alphas.end(), std::greater<double>());

      // Set interpolation knots and values
      _ipol.setQ2Values(q2s);
      _ipol.setAlphaSValues(alphas);
    }

    _calculated = true;
  }


}
