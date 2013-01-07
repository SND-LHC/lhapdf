#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/BicubicInterpolator.h"
#include <iostream>

namespace LHAPDF {


  double interpolateCubic(double T, double VL, double VDL, double VH, double VDH) {
    // Pre-calculate powers of T
    const double t2 = T*T;
    const double t3 = t2*T;

    // Calculate left point
    const double p0 = (2*t3 - 3*t2 + 1)*VL;
    const double m0 = (t3 - 2*t2 + T)*VDL;

    // Calculate right point
    const double p1 = (-2*t3 + 3*t2)*VH;
    const double m1 = (t3 - t2)*VDH;

    return p0 + m0 + p1 + m1;
  }


  // Provides d/dx at all grid locations
  /// @todo Convert to operate on a subgrid object instead
  double BicubicInterpolator::ddx(const double* xfs, size_t xidx, size_t q2idx) const {
    // Check for edge
    if (xidx == 0) {
      // Use forward difference
      return (xfs[pdf()->ptindex(xidx+1, q2idx)] - xfs[pdf()->ptindex(xidx, q2idx)]) / (pdf()->xKnots()[xidx+1]-pdf()->xKnots()[xidx]);
    } else if (xidx == pdf()->xKnots().size()-1) {
      // Use backward difference
      return (xfs[pdf()->ptindex(xidx, q2idx)] - xfs[pdf()->ptindex(xidx-1, q2idx)]) / (pdf()->xKnots()[xidx]-pdf()->xKnots()[xidx-1]);
    } else {
      // Use central difference
      double lddx = (xfs[pdf()->ptindex(xidx, q2idx)] - xfs[pdf()->ptindex(xidx-1, q2idx)]) / (pdf()->xKnots()[xidx]-pdf()->xKnots()[xidx-1]);
      double rddx = (xfs[pdf()->ptindex(xidx+1, q2idx)] - xfs[pdf()->ptindex(xidx, q2idx)]) / (pdf()->xKnots()[xidx+1]-pdf()->xKnots()[xidx]);
      return 0.5 * (lddx + rddx);
    }
  }


  double BicubicInterpolator::interpolateXQ2(int id, double x, double q2) const {
    /// Move this functionality to a helper which only operates on a subgrid, with interpolateXQ2 at a higher level to identify that subgrid

    // Index look-up
    size_t xidx = pdf()->xKnotLow(x);
    size_t q2idx = pdf()->q2KnotLow(q2);

    // Fractional parameters
    double tx = (x - pdf()->xKnots()[xidx]) / (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]);
    double tq = (q2 - pdf()->q2Knots()[q2idx]) / (pdf()->q2Knots()[q2idx+1] - pdf()->q2Knots()[q2idx]);

    // The xf value (sub)grid
    const double* xfs = pdf()->ptdata(id);

    // Points in Q2
    double vl = interpolateCubic(tx, xfs[pdf()->ptindex(xidx, q2idx)],
                                 ddx(xfs, xidx, q2idx) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                 xfs[pdf()->ptindex(xidx+1, q2idx)],
                                 ddx(xfs, xidx+1, q2idx) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );
    double vh = interpolateCubic(tx, xfs[pdf()->ptindex(xidx, q2idx+1)],
                                 ddx(xfs, xidx, q2idx+1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                 xfs[pdf()->ptindex(xidx+1, q2idx+1)],
                                 ddx(xfs, xidx+1, q2idx+1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );

    // Derivatives in Q2
    double vdl, vdh;

    if (q2idx == 0) {
      // Forward difference for lower q
      vdl = (vh - vl)/(pdf()->q2Knots()[q2idx+1]-pdf()->q2Knots()[q2idx]);
      // Central difference for higher q
      double vhh = interpolateCubic(tx, xfs[pdf()->ptindex(xidx, q2idx+2)],
                                    ddx(xfs, xidx, q2idx+2) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                    xfs[pdf()->ptindex(xidx+1, q2idx+2)],
                                    ddx(xfs, xidx+1, q2idx+2) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );
      vdh = 0.5 * vdl + 0.5 * (vhh - vh)/(pdf()->q2Knots()[q2idx+2] - pdf()->q2Knots()[q2idx+1]);
    }
    else if (q2idx+1 == pdf()->q2Knots().size()-1) {
      // Backward difference for higher q
      vdh = (vh - vl)/(pdf()->q2Knots()[q2idx+1]-pdf()->q2Knots()[q2idx]);
      // Central difference for lower q
      double vll = interpolateCubic(tx, xfs[pdf()->ptindex(xidx, q2idx-1)],
                                    ddx(xfs, xidx, q2idx-1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                    xfs[pdf()->ptindex(xidx+1, q2idx-1)],
                                    ddx(xfs, xidx+1, q2idx-1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );
      vdl = 0.5 * vdh + 0.5 * (vl - vll)/(pdf()->q2Knots()[q2idx] - pdf()->q2Knots()[q2idx-1]);
    }
    else {
      // Central difference for both q
      double vll = interpolateCubic( tx, xfs[pdf()->ptindex(xidx, q2idx-1)],
                                     ddx(xfs, xidx, q2idx-1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                     xfs[pdf()->ptindex(xidx+1, q2idx-1)],
                                     ddx(xfs, xidx+1, q2idx-1) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );
      vdl = 0.5 * (vh - vl)/(pdf()->q2Knots()[q2idx+1]-pdf()->q2Knots()[q2idx]) + 0.5 * (vl - vll)/(pdf()->q2Knots()[q2idx]-pdf()->q2Knots()[q2idx-1]);
      double vhh = interpolateCubic( tx, xfs[pdf()->ptindex(xidx, q2idx+2)],
                                     ddx(xfs, xidx, q2idx+2) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]),
                                     xfs[pdf()->ptindex(xidx+1, q2idx+2)],
                                     ddx(xfs, xidx+1, q2idx+2) * (pdf()->xKnots()[xidx+1] - pdf()->xKnots()[xidx]) );
      vdh = 0.5 * (vh - vl)/(pdf()->q2Knots()[q2idx+1] - pdf()->q2Knots()[q2idx]) + 0.5 * (vhh - vh)/(pdf()->q2Knots()[q2idx+2] - pdf()->q2Knots()[q2idx+1]);
    }

    vdl *= (pdf()->q2Knots()[q2idx+1]-pdf()->q2Knots()[q2idx]);
    vdh *= (pdf()->q2Knots()[q2idx+1]-pdf()->q2Knots()[q2idx]);

    double f_f = interpolateCubic( tq, vl, vdl, vh, vdh );

    return f_f;
  }


}
