#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Calculates first order derivative, dy/dt, as it appears in differential equation
  double AlphaS_ODE::_derivative(double t, double y, const std::vector<double>& beta ) {
    const double d0 = beta[0]/(-2*M_PI);
    const double d1 = beta[1]/(-4*M_PI*M_PI);
    const double d2 = beta[2]/(-64*M_PI*M_PI*M_PI);
    const double d = (d0*y*y) + (d1*y*y*y) + (d2*y*y*y*y);
    return d / t;
  }


  // Solves differential equation in alphaS using an implementation of RK4;
  double AlphaS_ODE::alphaS_Q2(double q2) {

    /// @todo Make these class members?
    // Step size
    double h = 0.001;
    // Fractional threshold to which we run in Q2 scale
    /// @todo Need a mechanism to shrink the steps if accuracy < stepsize?
    double accuracy = 0.001;

    // Check if we are going to run forward or backwards in energy scale towards target Q2.
    // Set h to step backwards if our target is less than where we are at.
    if (q2 < sqr(_mz)) h *= -1;

    // Number of active flavours used in beta function calculations
    const double nf = _nf_Q2(q2);
    const vector<double> bs = betas(nf);

    // Run in energy using RK4 algorithm until we are within our defined threshold energy
    double t = _mz; // starting point
    double y = _alphas_mz; // starting value
    double f1, f2, f3, f4, k1, k2, k3, k4;
    while (fabs(q2 - t)/q2 > accuracy) {
      // Increment based on the slope at the beginning of the interval using a simple Euler step
      f1 = _derivative(t, y, bs);
      k1 = h*f1;

      // Increments based on the slope at the midpoint of the interval
      f2 = _derivative(t + (h/2.0), y + (k1/2.0), bs);
      k2 = h*f2;

      f3 = _derivative(t + (h/2.0), y + (k2/2.0), bs);
      k3 = h*f3;

      // Increment based on the slope at the end of the interval
      f4 = _derivative(t + h, y + k3, bs);
      k4 = h*f4;

      // Calculate the weighted average of these increments and make the step
      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      t += h;
    }
    return y;
  }


}
