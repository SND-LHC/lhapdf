#include "LHAPDF/AlphaS.h"
#include <cmath>
#include <vector>
#include <iostream>

/// Masses of quarks in GeV.  Used to calculate the number of quarks that are active at a given energy range Q2
/// @todo Put into set metadata
static const double quarkMasses[] = {0.0017, 0.0041, 0.1, 1.29, 4.1, 172};
/// @todo Allow for general lambdaQCD, set to 0.4 for now
static double lambdaQCD = 0.4;


/// Calculate the number of active quark flavours at energy scale Q2
int AlphaS::getNfAtQ2(double q2) {
  int nf = 2;
  for ( ; nf < 6; ++nf) {
    if (q2 < sqr(quarkMasses[nf])) break; // check against *next* mass
  }
  return nf; // max value = 5
}


/// Calculate beta functions given the number of active flavours
std::vector<double> AlphaS::getBetas(int nf) {
  std::vector<double> v;
  double beta0 = 11 - (2 / 3.0)*nf;
  double beta1 = 51 - (19 / 3.0)*nf;
  double beta2 = 2857 - (5033 / 9.0)*nf + (325 / 27.0)*nf*nf;
  v.push_back(beta0);
  v.push_back(beta1);
  v.push_back(beta2);
  return v;
}


/// @todo WTF?
double AlphaS::meta(double q2) {
  return 0.0;
}


/// Calculate alpha_s(Q2) by an analytic approximation
double AlphaS::analytic(double q2, double lambdaQCD, int order=2) {
  // Get beta coeffs for the number of active (above threshold) quark flavours at energy Q
  const int nf = getNfAtQ2( Q2 );
  const std::vector<double> beta = getBetas(nf);
  const double beta02 = beta[0] * beta[0];
  const double beta12 = beta[1] * beta[1];

  /// @todo Use both 4 and 5 flavour LambdaQCD definitions depending on m_b threshold

  // Pre-calculate ln(Q2/lambdaQCD) and expansion term y = 1/ln(Q2/lambdaQCD)
  const double x = q2 / lambdaQCD;
  const double lnx = log(x);
  const double lnlnx = log(lnx);
  const double y = 1 / lnx;

  // Calculate expansion terms Ay[a_0*y^0 - a_1*y^1 + By^2[a_20 + a_21 - a_22]] where a_ij are series coefficients
  const double A = 4 * M_PI / beta[0];
  const double a_0 = 1;
  const double a_1 = - 2 * beta[1] * lnlnx / beta02;
  const double B = 4 * beta12 / (beta02 * beta02);
  const double a_20 = (lnlnx - 0.5) * (lnlnx - 0.5);
  const double a_21 = beta[2] * beta[0] / (8 * beta12);
  const double a_22 = 5 / 4.0;

  // Combine terms
  double tmp = a_0;
  if (order > 0) tmp -= a_1 * y;
  if (order > 1) tmp += B * y*y * (a_20 + a_21 - a_22);
  const double alphaS = A * y * tmp;
  return alphaS;
}


// Calculates first order derivative, dy/dt, as it appears in differential equation
double AlphaS::derivative( double t, double y, std::vector<double> beta ) {
  const double d0 = beta[0]/(-2*M_PI);
  const double d1 = beta[1]/(-4*M_PI*M_PI);
  const double d2 = beta[2]/(-64*M_PI*M_PI*M_PI);
  const double derivative = (d0*y*y) + (d1*y*y*y) + (d2*y*y*y*y);
  derivative /= t;
  return derivative;
}


// Solves differential equation in alphaS using an implementation of RK4;
double AlphaS::numerical(double q2) {

  /// @todo Make these class members
  // Step size
  double h = 0.001;
  // Fractional threshold to which we run in energy scale
  double accuracy = 0.1;
  // Mass of the Z-boson
  double m_z = 91.1876;
  // Value of alpha_s(MZ)
  double alphaSmZ = 0.118;


  // Check if we are going to run forward or backwards in energy scale towards target Q2
  // Set h to step backwards if our target is less than where we are at
  if (q2 < sqr(m_z)) h *= -1;

  // Number of active flavours used in beta function calculations
  const double nf = getNfAtQ2(q2);
  const std::vector<double> beta = getBetas(nf);

  // Run in energy using RK4 algorithm until we are within our defined threshold energy
  double t = m_z; // starting point
  double y = alphaSmZ; // starting value
  double f1, f2, f3, f4, k1, k2, k3, k4;
  while (fabs(q2 - t)/q2 > accuracy) {
    // Increment based on the slope at the beginning of the interval using a simple Euler step
    f1 = derivative( t, y, beta );
    k1 = h*f1;

    // Increments based on the slope at the midpoint of the interval
    f2 = derivative( t + (h/2.0), y + (k1/2.0), beta );
    k2 = h*f2;

    f3 = derivative( t + (h/2.0), y + (k2/2.0), beta );
    k3 = h*f3;

    // Increment based on the slope at the end of the interval
    f4 = derivative( t + h, y + k3, beta );
    k4 = h*f4;

    // Calculate the weighted average of these increments and make the step
    y += (1/6.0)*(k1 + (2*k2) + (2*k3) + k4);
    t += h;
  }
  return y;
}
