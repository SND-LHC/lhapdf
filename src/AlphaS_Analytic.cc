#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  /// Calculate alpha_s(Q2) by an analytic approximation
  double AlphaS_Analytic::alphaS_Q2(double q2) {
    // Get beta coeffs for the number of active (above threshold) quark flavours at energy Q
    const int nf = _nf_Q2(q2);
    const std::vector<double> beta = betas(nf);
    const double beta02 = beta[0] * beta[0];
    const double beta12 = beta[1] * beta[1];

    /// @todo Use both 4 and 5 flavour LambdaQCD definitions depending on m_b threshold
    const double lambdaQCD = _lambda4;

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
    if (_order > 0) tmp -= a_1 * y;
    if (_order > 1) tmp += B * y*y * (a_20 + a_21 - a_22);
    const double alphaS = A * y * tmp;
    return alphaS;
  }


}
