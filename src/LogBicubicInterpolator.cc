// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/LogBicubicInterpolator.h"
#include <iostream>

namespace LHAPDF {


  namespace { // Unnamed namespace

    /// One-dimensional linear interpolation for y(x)
    inline double _interpolateLinear(double x, double xl, double xh, double yl, double yh)	{
      assert(x >= xl);
      assert(xh >= x);
      return yl + (x - xl) / (xh - xl) * (yh - yl);
    }

    /// One-dimensional cubic interpolation
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


    /// Calculate adjacent d(xf)/dx at all grid locations for fixed iq2
    ///
    /// @todo Store pre-cached dlogxs, dlogq2s on subgrids, to replace these denominators? Any real speed gain for the extra memory?
    double _dxf_dlogx(const KnotArray1F& subgrid, size_t ix, size_t iq2) {
      const size_t nxknots = subgrid.xs().size();
      if (ix != 0 && ix != nxknots-1) { //< If central, use the central difference
        /// @note We evaluate the most likely condition first to help compiler branch prediction
        const double lddx = (subgrid.xf(ix, iq2) - subgrid.xf(ix-1, iq2)) / (subgrid.logxs()[ix] - subgrid.logxs()[ix-1]);
        const double rddx = (subgrid.xf(ix+1, iq2) - subgrid.xf(ix, iq2)) / (subgrid.logxs()[ix+1] - subgrid.logxs()[ix]);
        return (lddx + rddx) / 2.0;
      } else if (ix == 0) { //< If at leftmost edge, use forward difference
        return (subgrid.xf(ix+1, iq2) - subgrid.xf(ix, iq2)) / (subgrid.logxs()[ix+1] - subgrid.logxs()[ix]);
      } else if (ix == nxknots-1) { //< If at rightmost edge, use backward difference
        return (subgrid.xf(ix, iq2) - subgrid.xf(ix-1, iq2)) / (subgrid.logxs()[ix] - subgrid.logxs()[ix-1]);
      } else {
        throw LogicError("We shouldn't be able to get here!");
      }
    }

  }


  // Static method for cache acquisition
  LogBicubicInterpolator::XQ2Cache& LogBicubicInterpolator::_getCache(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) {
    /// Meyers Singleton is automatically thread-local since block scope:
    /// https://www.modernescpp.com/index.php/thread-safe-initialization-of-a-singleton
    static XQ2Cache instance;

    /// @todo Multi-cache as maps on the grid hashes!
    /// @todo Fuzzy testing on x and Q2?

    // Update the x cache
    const size_t xhash = subgrid.xhash();
    const bool xok = instance.x == x;
    cout << instance.xhash << " " << xhash << endl;
    const bool ixok = instance.ix == ix && instance.xhash == xhash;
    if (!xok) {
      instance.logx = log(x);
    }
    if (!ixok) {
      instance.xhash = xhash;
      instance.dlogx_1 = subgrid.logxs()[ix+1] - subgrid.logxs()[ix];
    }
    if (!xok || !ixok) {
      instance.tlogx = (instance.logx - subgrid.logxs()[ix]) / instance.dlogx_1;
    }

    // Update the Q2 cache
    const size_t q2hash = subgrid.q2hash();
    const bool q2ok = instance.q2 == q2;
    const bool iq2ok = instance.iq2 == iq2 && instance.q2hash == q2hash;
    if (!q2ok) {
      instance.logq2 = log(q2);
    }
    if (!iq2ok) {
      instance.q2hash = q2hash;
      instance.dlogq_0 = (iq2 != 0) ? subgrid.logq2s()[iq2] - subgrid.logq2s()[iq2-1] : -1; //< Don't evaluate (or use) if iq2-1 < 0
      instance.dlogq_1 = subgrid.logq2s()[iq2+1] - subgrid.logq2s()[iq2];
      instance.dlogq_2 = (iq2+2 != subgrid.xsize()) ? subgrid.logq2s()[iq2+2] - subgrid.logq2s()[iq2+1] : -1; //< Don't evaluate (or use) if iq2+2 > iq2max
      // instance.dlogq_2 = (iq2+1 != iq2max) ? subgrid.logq2s()[iq2+2] - subgrid.logq2s()[iq2+1] : -1; //< Don't evaluate (or use) if iq2+2 > iq2max
    }
    if (!q2ok || !iq2ok) {
      instance.tlogq = (instance.logq2 - subgrid.logq2s()[iq2]) / instance.dlogq_1;
    }

    return instance;
  }


  double LogBicubicInterpolator::_interpolateXQ2(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) const {
    // Raise an error if there are too few knots even for a linear fall-back
    const size_t nxknots = subgrid.xsize();
    const size_t nq2knots = subgrid.q2size();
    if (nxknots < 4)
      throw GridError("PDF subgrids are required to have at least 4 x-knots for use with LogBicubicInterpolator");
    if (nq2knots < 2)
      throw GridError("PDF subgrids are required to have at least 2 Q-knots for use with LogBicubicInterpolator");

    // Check x and q index ranges -- we always need i and i+1 indices to be valid
    const size_t ixmax = nxknots - 1;
    const size_t iq2max = nq2knots - 1;
    if (ix+1 > ixmax) // also true if ix is off the end
      throw GridError("Attempting to access an x-knot index past the end of the array, in linear fallback mode");
    if (iq2+1 > iq2max) // also true if iq2 is off the end
      throw GridError("Attempting to access an Q-knot index past the end of the array, in linear fallback mode");

    // Check, update, and get the cache
    const XQ2Cache& cache = _getCache(subgrid, x, ix, q2, iq2);
    const double logx = cache.logx;
    const double logq2 = cache.logq2;

    // Fall back to LogBilinearInterpolator if either 2 or 3 Q-knots
    if (nq2knots < 4) {
      // First interpolate in x
      const double logx0 = subgrid.logxs()[ix];
      const double logx1 = subgrid.logxs()[ix+1];
      const double f_ql = _interpolateLinear(logx, logx0, logx1, subgrid.xf(ix, iq2), subgrid.xf(ix+1, iq2));
      const double f_qh = _interpolateLinear(logx, logx0, logx1, subgrid.xf(ix, iq2+1), subgrid.xf(ix+1, iq2+1));
      // Then interpolate in Q2, using the x-ipol results as anchor points
      return _interpolateLinear(logq2, subgrid.logq2s()[iq2], subgrid.logq2s()[iq2+1], f_ql, f_qh);
    }
    // else proceed with cubic interpolation:

    // Pre-calculate parameters
    /// @todo Remove this redundant step: it's now just making aliases for the cache variables
    const double& dlogx_1 = cache.dlogx_1;
    const double& tlogx = cache.tlogx;
    const double& dlogq_0 = cache.dlogq_0; //< Don't evaluate (or use) if iq2-1 < 0
    const double& dlogq_1 = cache.dlogq_1;
    const double& dlogq_2 = cache.dlogq_2; //< Don't evaluate (or use) if iq2+2 > iq2max
    const double& tlogq = cache.tlogq;

    /// @todo Statically pre-compute the whole nx * nq gradiant array? I.e. _dxf_dlogx for all points in all subgrids. Memory ~doubling :-/ Could cache them as they are used...

    // Points in Q2
    double vl = _interpolateCubic(tlogx, subgrid.xf(ix, iq2), _dxf_dlogx(subgrid, ix, iq2) * dlogx_1,
                                         subgrid.xf(ix+1, iq2), _dxf_dlogx(subgrid, ix+1, iq2) * dlogx_1);
    double vh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+1), _dxf_dlogx(subgrid, ix, iq2+1) * dlogx_1,
                                         subgrid.xf(ix+1, iq2+1), _dxf_dlogx(subgrid, ix+1, iq2+1) * dlogx_1);

    // Derivatives in Q2
    double vdl, vdh;
    if (iq2 > 0 && iq2+1 < iq2max) {
      // Central difference for both q
      /// @note We evaluate the most likely condition first to help compiler branch prediction
      double vll = _interpolateCubic(tlogx, subgrid.xf(ix, iq2-1), _dxf_dlogx(subgrid, ix, iq2-1) * dlogx_1,
                                            subgrid.xf(ix+1, iq2-1), _dxf_dlogx(subgrid, ix+1, iq2-1) * dlogx_1);
      vdl = ( (vh - vl)/dlogq_1 + (vl - vll)/dlogq_0 ) / 2.0;
      double vhh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+2), _dxf_dlogx(subgrid, ix, iq2+2) * dlogx_1,
                                            subgrid.xf(ix+1, iq2+2), _dxf_dlogx(subgrid, ix+1, iq2+2) * dlogx_1);
      vdh = ( (vh - vl)/dlogq_1 + (vhh - vh)/dlogq_2 ) / 2.0;
    }
    else if (iq2 == 0) {
      // Forward difference for lower q
      vdl = (vh - vl) / dlogq_1;
      // Central difference for higher q
      double vhh = _interpolateCubic(tlogx, subgrid.xf(ix, iq2+2), _dxf_dlogx(subgrid, ix, iq2+2) * dlogx_1,
                                            subgrid.xf(ix+1, iq2+2), _dxf_dlogx(subgrid, ix+1, iq2+2) * dlogx_1);
      vdh = (vdl + (vhh - vh)/dlogq_2) / 2.0;
    }
    else if (iq2+1 == iq2max) {
      // Backward difference for higher q
      vdh = (vh - vl) / dlogq_1;
      // Central difference for lower q
      double vll = _interpolateCubic(tlogx, subgrid.xf(ix, iq2-1), _dxf_dlogx(subgrid, ix, iq2-1) * dlogx_1,
                                            subgrid.xf(ix+1, iq2-1), _dxf_dlogx(subgrid, ix+1, iq2-1) * dlogx_1);
      vdl = (vdh + (vl - vll)/dlogq_0) / 2.0;
    }
    else throw LogicError("We shouldn't be able to get here!");

    vdl *= dlogq_1;
    vdh *= dlogq_1;
    return _interpolateCubic(tlogq, vl, vdl, vh, vdh);
  }


}
