// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/LogBicubicInterpolator.h"
#include <iostream>
#include <mutex>
#include <thread>

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


  // Static method for x cache acquisition
  LogBicubicInterpolator::XCache& LogBicubicInterpolator::_getCacheX(const KnotArray1F& subgrid, double x, size_t ix) {
    // static mutex xmutex;
    static map<thread::id,XCachesMap> xcachemaps; //< thread-safe Meyers Singleton

    // Get the thread-local caches
    const thread::id tid = this_thread::get_id();
    XCacheMap& xcaches = xcachemaps[tid]; ///< @todo Need a lock for initialisation?

    // Get the subgrid-specific x cache
    const size_t xhash = subgrid.xhash();
    XCache& xcache = xcaches[xhash];

    // Check the multi-level cache
    /// @todo Make cache multi-level
    /// @todo Cache more ipol-weight variables?
    /// @todo Fuzzy testing on x?
    /// @todo Push back to an earlier stage, to avoid-refinding ix for same x and same grid
    for (size_t i = 0; i < xcache.N; ++i) {
      const size_t j = (xcache.ilast - i) % xcache.N; //< step backwards, deeper into history
      const bool xok = xcache.x[j] == x;
      // const bool ixok = xcache.ix == ix;
      assert(xcache.ix[j] == ix); //< guaranteed to be same subgrid, so ix *must* be the same
      if (!xok) continue;


    // Update the x cache
    if (!xok) {
      if (!xok) {
        xcache.logx = log(x);
      }
      if (!ixok) {
        // xcache.xhash = xhash;
        xcache.dlogx_1 = subgrid.logxs()[ix+1] - subgrid.logxs()[ix];
      }
      xcache.tlogx = (xcache.logx - subgrid.logxs()[ix]) / xcache.dlogx_1;
    }

    return xcache;
  }


  // Static method for Q2 cache acquisition
  LogBicubicInterpolator::Q2Cache& LogBicubicInterpolator::_getCacheQ2(const KnotArray1F& subgrid, double q2, size_t iq2) {
    // static mutex q2mutex;
    static map<thread::id,Q2CacheMap> q2cachemaps; //< thread-safe Meyers Singleton

    // Get the thread-local caches
    const thread::id tid = this_thread::get_id();
    Q2CacheMap& q2caches = q2cachemaps[tid]; ///< @todo Need a lock for initialisation?

    // Get and check the subgrid-specific Q2 cache
    /// @todo Make cache multi-level
    /// @todo Fuzzy testing on Q2?
    /// @todo Cache more ipol-weight variables?
    const size_t q2hash = subgrid.q2hash();
    Q2Cache& q2cache = q2caches[q2hash];
    const bool q2ok = q2cache.q2 == q2;
    const bool iq2ok = q2cache.iq2 == iq2; // && q2cache.q2hash == q2hash;

    // Update the Q2 cache
    if (!q2ok || !iq2ok) {
      if (!q2ok) {
        q2cache.logq2 = log(q2);
      }
      if (!iq2ok) {
        q2cache.q2hash = q2hash;
        q2cache.dlogq_0 = (iq2 != 0) ? subgrid.logq2s()[iq2] - subgrid.logq2s()[iq2-1] : -1; //< Don't evaluate (or use) if iq2-1 < 0
        q2cache.dlogq_1 = subgrid.logq2s()[iq2+1] - subgrid.logq2s()[iq2];
        q2cache.dlogq_2 = (iq2+2 != subgrid.xsize()) ? subgrid.logq2s()[iq2+2] - subgrid.logq2s()[iq2+1] : -1; //< Don't evaluate (or use) if iq2+2 > iq2max
      }
      q2cache.tlogq = (q2cache.logq2 - subgrid.logq2s()[iq2]) / q2cache.dlogq_1;
    }

    return q2cache;
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
    const XCache& xcache = _getCacheX(subgrid, x, ix);
    const Q2Cache& q2cache = _getCacheQ2(subgrid, q2, iq2);
    const double logx = xcache.logx;
    const double logq2 = q2cache.logq2;

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
    const double& dlogx_1 = xcache.dlogx_1;
    const double& tlogx = xcache.tlogx;
    const double& dlogq_0 = q2cache.dlogq_0; //< Don't evaluate (or use) if iq2-1 < 0
    const double& dlogq_1 = q2cache.dlogq_1;
    const double& dlogq_2 = q2cache.dlogq_2; //< Don't evaluate (or use) if iq2+2 > iq2max
    const double& tlogq = q2cache.tlogq;

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
