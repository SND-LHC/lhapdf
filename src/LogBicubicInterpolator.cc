#include "LHAPDF/LogBicubicInterpolator.h"
#include <iostream>

namespace LHAPDF {


  namespace { // Unnamed namespace

    // One-dimensional cubic interpolation
    inline double _interpolateCubic(double T, double VL, double VDL, double VH, double VDH) {
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
    double _ddlogx(const KnotArray1F& subgrid, size_t ix, size_t iq2) {
      /// @todo Re-order this if so that branch prediction will favour the "normal" central case
      if (ix == 0) { //< If at leftmost edge, use forward difference
        return (subgrid.xf(ix+1, iq2) - subgrid.xf(ix, iq2)) / (log(subgrid.xs()[ix+1]) - log(subgrid.xs()[ix]));
      } else if (ix == subgrid.xs().size() - 1) { //< If at rightmost edge, use backward difference
        return (subgrid.xf(ix, iq2) - subgrid.xf(ix-1, iq2)) / (log(subgrid.xs()[ix]) - log(subgrid.xs()[ix-1]));
      } else { //< If central, use the central difference
        const double lddx = (subgrid.xf(ix, iq2) - subgrid.xf(ix-1, iq2)) / (log(subgrid.xs()[ix]) - log(subgrid.xs()[ix-1]));
        const double rddx = (subgrid.xf(ix+1, iq2) - subgrid.xf(ix, iq2)) / (log(subgrid.xs()[ix+1]) - log(subgrid.xs()[ix]));
        return (lddx + rddx) / 2.0;
      }
    }

  }



  double LogBicubicInterpolator::_interpolateXQ2(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) const {
    /// @todo Allow interpolation right up to the borders of the grid in Q2 and x... the last inter-knot range is currently broken

    /// @todo Also treat the x top/bottom edges carefully, cf. the Q2 ones

    const double logx = log(x);
    const double logq2 = log(q2);

    // Distance parameters
    const double dlogx = log(subgrid.xs()[ix+1]) - log(subgrid.xs()[ix]);
    const double tlogx = (logx - log(subgrid.xs()[ix])) / dlogx;
    /// @todo Only compute these if the +1 and +2 indices are guaranteed to be valid
    const double dlogq_0 = log(subgrid.q2s()[iq2]) - log(subgrid.q2s()[iq2-1]);
    const double dlogq_1 = log(subgrid.q2s()[iq2+1]) - log(subgrid.q2s()[iq2]);
    const double dlogq_2 = log(subgrid.q2s()[iq2+2]) - log(subgrid.q2s()[iq2+1]);
    const double dlogq = dlogq_1;
    const double tlogq = (logq2 - log(subgrid.q2s()[iq2])) / dlogq;

    // Points in Q2
    double vl = _interpolateCubic(tlogx, subgrid.xf(ix, iq2), _ddlogx(subgrid, ix, iq2) * dlogx,
                                         subgrid.xf(ix+1, iq2), _ddlogx(subgrid, ix+1, iq2) * dlogx);
    double vh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+1), _ddlogx(subgrid, ix, iq2+1) * dlogx,
                                         subgrid.xf(ix+1, iq2+1), _ddlogx(subgrid, ix+1, iq2+1) * dlogx);

    // Derivatives in Q2
    double vdl, vdh;
    if (iq2 == 0) {
      // Forward difference for lower q
      vdl = (vh - vl) / dlogq_1;
      // Central difference for higher q
      double vhh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+2), _ddlogx(subgrid, ix, iq2+2) * dlogx,
                                            subgrid.xf(ix+1, iq2+2), _ddlogx(subgrid, ix+1, iq2+2) * dlogx);
      vdh = (vdl + (vhh - vh)/dlogq_2) / 2.0;
    }
    else if (iq2+1 == subgrid.q2s().size()-1) {
      // Backward difference for higher q
      vdh = (vh - vl) / dlogq_1;
      // Central difference for lower q
      double vll = _interpolateCubic(tlogx, subgrid.xf(ix, iq2-1), _ddlogx(subgrid, ix, iq2-1) * dlogx,
                                            subgrid.xf(ix+1, iq2-1), _ddlogx(subgrid, ix+1, iq2-1) * dlogx);
      vdl = (vdh + (vl - vll)/dlogq_0) / 2.0;
    }
    else {
      // Central difference for both q
      double vll = _interpolateCubic(tlogx, subgrid.xf(ix, iq2-1), _ddlogx(subgrid, ix, iq2-1) * dlogx,
                                            subgrid.xf(ix+1, iq2-1), _ddlogx(subgrid, ix+1, iq2-1) * dlogx);
      vdl = ( (vh - vl)/dlogq_1 + (vl - vll)/dlogq_0 ) / 2.0;
      double vhh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+2), _ddlogx(subgrid, ix, iq2+2) * dlogx,
                                            subgrid.xf(ix+1, iq2+2), _ddlogx(subgrid, ix+1, iq2+2) * dlogx);
      vdh = ( (vh - vl)/dlogq_1 + (vhh - vh)/dlogq_2 ) / 2.0;
    }

    vdl *= dlogq;
    vdh *= dlogq;

    return _interpolateCubic(tlogq, vl, vdl, vh, vdh);
  }


}
