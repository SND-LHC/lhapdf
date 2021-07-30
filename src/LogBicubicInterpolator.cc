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

  }
  
  double LogBicubicInterpolator::_interpolateXQ2(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, int id) const {
    
    // Raise an error if there are too few knots even for a linear fall-back
    const size_t nxknots = grid.xsize();
    const size_t nq2knots = grid.q2size();

    // MK: do you really need different number of knots in the directions?
    //   probably should be <2 for both methods, and fall back to linear in both cases.
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

    // Fall back to LogBilinearInterpolator if either 2 or 3 Q-knots
    // these two lines make for roughly 20% of the computing time
    const double logx  = log(x);
    const double logq2 = log(q2);

    if (nq2knots < 4) {
      // First interpolate in x
      const double logx0 = grid.logxs(ix);
      const double logx1 = grid.logxs(ix+1);
      const double f_ql = _interpolateLinear(logx, logx0, logx1, grid.xf(ix, iq2, id),   grid.xf(ix+1, iq2, id));
      const double f_qh = _interpolateLinear(logx, logx0, logx1, grid.xf(ix, iq2+1, id), grid.xf(ix+1, iq2+1, id));
      // Then interpolate in Q2, using the x-ipol results as anchor points
      return _interpolateLinear(logq2, grid.logq2s(iq2), grid.logq2s(iq2+1), f_ql, f_qh);
    }
    // else proceed with cubic interpolation:

    // Pre-calculate parameters
    const double& dlogx_1 = grid.logxs(ix+1) - grid.logxs(ix);
    const double& tlogx   = (logx - grid.logxs(ix)) / dlogx_1;    
    const double& dlogq_0 = (iq2 != 0) ? grid.logq2s(iq2) - grid.logq2s(iq2-1) : -1; //< Don't evaluate (or use) if iq2-1 < 0
    const double& dlogq_1 = grid.logq2s(iq2+1) - grid.logq2s(iq2);    
    const double& dlogq_2 = (iq2+2 != grid.q2size()) ? grid.logq2s(iq2+2) - grid.logq2s(iq2+1) : -1; //< Don't evaluate (or use) if iq2+2 > iq2max
    const double& tlogq   = (logq2 - grid.logq2s(iq2)) / dlogq_1;

    // Points in Q2
    double vl = _interpolateCubic(tlogx, grid.xf(ix, iq2, id), grid.dxf(ix, iq2, id) * dlogx_1,
				  grid.xf(ix+1, iq2, id), grid.dxf(ix+1, iq2, id) * dlogx_1);
    double vh = _interpolateCubic(tlogx, grid.xf(ix, iq2+1, id), grid.dxf(ix, iq2+1, id) * dlogx_1,
				  grid.xf(ix+1, iq2+1, id), grid.dxf(ix+1, iq2+1, id) * dlogx_1);


    const bool q2_lower = ( (iq2 == 0) || (grid.q2s(iq2) == grid.q2s(iq2-1)));
    const bool q2_upper = ( (iq2 == grid.q2size() -1) || (grid.q2s(iq2+1) == grid.q2s(iq2+2)) );

    // Derivatives in Q2
    double vdl, vdh;
    if (!q2_lower && !q2_upper){
      // Central difference for both q
      /// @note We evaluate the most likely condition first to help compiler branch prediction
      double vll = _interpolateCubic(tlogx, grid.xf(ix, iq2-1, id), grid.dxf(ix, iq2-1, id) * dlogx_1,
				     grid.xf(ix+1, iq2-1, id), grid.dxf(ix+1, iq2-1, id) * dlogx_1);
      vdl = ( (vh - vl)/dlogq_1 + (vl - vll)/dlogq_0 ) / 2.0;
      double vhh = _interpolateCubic(tlogx, grid.xf(ix, iq2+2, id), grid.dxf(ix, iq2+2, id) * dlogx_1,
				     grid.xf(ix+1, iq2+2, id), grid.dxf(ix+1, iq2+2, id) * dlogx_1);
      vdh = ( (vh - vl)/dlogq_1 + (vhh - vh)/dlogq_2 ) / 2.0;
    }
    else if (q2_lower) {
      // Forward difference for lower q
      vdl = (vh - vl) / dlogq_1;
      // Central difference for higher q
      double vhh = _interpolateCubic(tlogx, grid.xf(ix, iq2+2, id), grid.dxf(ix, iq2+2, id) * dlogx_1,
				     grid.xf(ix+1, iq2+2, id), grid.dxf(ix+1, iq2+2, id) * dlogx_1);
      vdh = (vdl + (vhh - vh)/dlogq_2) / 2.0;
    }
    else if (q2_upper) {
      // Backward difference for higher q
      vdh = (vh - vl) / dlogq_1;
      // Central difference for lower q
      double vll = _interpolateCubic(tlogx, grid.xf(ix, iq2-1, id), grid.dxf(ix, iq2-1, id) * dlogx_1,
				     grid.xf(ix+1, iq2-1, id), grid.dxf(ix+1, iq2-1, id) * dlogx_1);
      vdl = (vdh + (vl - vll)/dlogq_0) / 2.0;
    }
    else throw LogicError("We shouldn't be able to get here!");

    vdl *= dlogq_1;
    vdh *= dlogq_1;
    return _interpolateCubic(tlogq, vl, vdl, vh, vdh);
  }


}
