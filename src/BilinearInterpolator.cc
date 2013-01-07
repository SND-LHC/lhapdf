#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/BilinearInterpolator.h"
#include <algorithm>
#include <cassert>

using namespace std;

namespace LHAPDF {


  /// One-dimensional linear interpolation for y(x)
  inline double interpolateLinear(double x, double xl, double xh, double yl, double yh)	{
    assert(x > xl);
    assert(xh > x);
    return yl + (x - xl) / (xh - xl) * (yh - yl);
  }


  double BilinearInterpolator::interpolateXQ2(int id, double x, double q2) const {
    // Calculate interpolation along X axis
    const vector<double>& xpoints = pdf()->xKnots();
    const vector<double>& qpoints = pdf()->q2Knots();
    size_t ix = pdf()->xKnotLow(x);
    size_t iq2 = pdf()->q2KnotLow(q2);

    /// @todo Cache this lookup... clean?

    // First interpolate in x
    const double* f = pdf()->ptdata(id);
    const double f_ql = interpolateLinear(x, xpoints[ix], xpoints[ix+1],
                                          f[pdf()->ptindex(ix, iq2)], f[pdf()->ptindex(ix+1, iq2)]);
    const double f_qh = interpolateLinear(x, xpoints[ix], xpoints[ix+1],
                                          f[pdf()->ptindex(ix, iq2+1)], f[pdf()->ptindex(ix+1, iq2+1)]);

    // Then interpolate in Q2, using the x-ipol results as anchor points
    return interpolateLinear(q2, qpoints[iq2], qpoints[iq2+1], f_ql, f_qh);
  }


}
