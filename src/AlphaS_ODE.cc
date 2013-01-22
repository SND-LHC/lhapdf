#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  namespace { // unnamed namespace

    // Calculate first order derivative, dy/dt, as it appears in the differential equation
    double _derivative(double t, double y, const std::vector<double>& beta ) {
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
    double h = 0.01;
    // Fractional threshold to which we run in Q2 scale
    /// @todo Need a mechanism to shrink the steps if accuracy < stepsize?
    /// @todo Is the accuracy in Q or Q2? Q is better, I think.
    const double faccuracy = 0.001;

    // Calculated absolute tolerance around the target Q2.
    /// @todo Would specifying the accuracy in abs rather than rel terms be better?
    const double accuracy = faccuracy * q2;

    // Number of active flavours used in beta function calculations
    /// @todo This needs to happen *inside* the running loop
    const int nf = nf_Q2(q2);
    const vector<double> bs = betas(nf);

    /// @todo Running (i.e. t and h) is in Q rather than Q2? Check

    /// @todo Use caching of the last target so we don't always have to start the running from MZ?
    ///       Then choose to start from last or MZ, whichever is closer.

    // Run in energy using RK4 algorithm until we are within our defined threshold energy
    double t = _mz; // starting point
    double y = _alphas_mz; // starting value
    /// @todo Units of q2 and t are different. Which quantity are we evolving in?
    while (fabs(q2 - t) > accuracy) {
      /// Mechanism to shrink the steps if accuracy < stepsize and we are close to Q2
      if (fabs(h) > accuracy && fabs(q2 - t)/h < 10) h = accuracy/2.0;

      // Check if we are going to run forward or backwards in energy scale towards target Q2.
      /// @todo C++11's std::copysign would be perfect here
      if ((q2 < t && sgn(h) > 0) || (q2 > t && sgn(h) < 0)) h *= -1;

      // cout << t << " -> " << q2 << ", in steps of size " << h << endl;

      // Increment based on the slope at the beginning of the interval using a simple Euler step
      const double k1 = h * _derivative(t, y, bs);

      // Increments based on the slope at the midpoint of the interval
      const double k2 = h * _derivative(t + h/2.0, y + k1/2.0, bs);

      const double k3 = h * _derivative(t + h/2.0, y + k2/2.0, bs);

      // Increment based on the slope at the end of the interval
      const double k4 = h * _derivative(t + h, y + k3, bs);

      // Calculate the weighted average of these increments and make the step
      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      t += h;
    }
    return y;
  }


}
