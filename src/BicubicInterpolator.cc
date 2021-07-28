// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/BicubicInterpolator.h"
#include <iostream>

namespace LHAPDF {


  namespace { // Unnamed namespace

    // One-dimensional linear interpolation for y(x)
    inline double _interpolateLinear(double x, double xl, double xh, double yl, double yh)	{
      assert(x >= xl);
      assert(xh >= x);
      return yl + (x - xl) / (xh - xl) * (yh - yl);
    }

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
    double _ddx(const KnotArray& grid, size_t ix, size_t iq2, int id) {
      /// @todo Re-order this if so that branch prediction will favour the "normal" central case
      if (ix == 0) { //< If at leftmost edge, use forward difference
        return (grid.xf(ix+1, iq2, id) - grid.xf(ix, iq2, id)) / (grid.xs(ix+1) - grid.xs(ix));
      } else if (ix == grid.xsize() - 1) { //< If at rightmost edge, use backward difference
        return (grid.xf(ix, iq2, id) - grid.xf(ix-1, iq2, id)) / (grid.xs(ix) - grid.xs(ix-1));
      } else { //< If central, use the central difference
        const double lddx = (grid.xf(ix, iq2, id) - grid.xf(ix-1, iq2, id)) / (grid.xs(ix) - grid.xs(ix-1));
        const double rddx = (grid.xf(ix+1, iq2, id) - grid.xf(ix, iq2, id)) / (grid.xs(ix+1) - grid.xs(ix));
        return (lddx + rddx) / 2.0;
      }
    }

  }


  
  double BicubicInterpolator::_interpolateXQ2(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, int id) const {

    // MK: if they have at least 2 knots, falls back to linear interpolator
    if (grid.xsize() < 4)
      throw GridError("PDF subgrids are required to have at least 4 x-knots for use with BicubicInterpolator");    
    if (grid.q2size() < 4) 
      throw GridError("PDF subgrids are required to have at least 4 Q2-knots for use with BicubicInterpolator");
    
    /// @todo Allow interpolation right up to the borders of the grid in Q2 and x... the last inter-knot range is currently broken
    /// @todo Also treat the x top/bottom edges carefully, cf. the Q2 ones

    // check edges, including internal discontinuity
    const bool q2_lower = ( (iq2 == 0) || (grid.q2s(iq2) == grid.q2s(iq2-1)));
    const bool q2_upper = ( (iq2 == grid.q2size() -1) || (grid.q2s(iq2+1) == grid.q2s(iq2+2)) );
    //const bool ix_lower = ( (ix == 0) );
    //const bool ix_upper = ( (ix == grid.xsize()) );
    
    // Distance parameters
    const double dx = grid.xs(ix+1) - grid.xs(ix);
    const double tx = (x - grid.xs(ix)) / dx;
    /// @todo Only compute these if the +1 and +2 indices are guaranteed to be valid
    // i.e. check if that is in range, and there is no discontinuitie there
    const double dq_0 = grid.q2s(iq2  ) - grid.q2s(iq2-1);
    const double dq_1 = grid.q2s(iq2+1) - grid.q2s(iq2  );
    const double dq_2 = grid.q2s(iq2+2) - grid.q2s(iq2+1);
    const double dq = dq_1;
    const double tq = (q2 - grid.q2s(iq2)) / dq;

    // Points in Q2
    double vl = _interpolateCubic(tx, grid.xf(ix,   iq2,   id), _ddx(grid, ix, iq2, id) * dx,
				  grid.xf(ix+1, iq2,   id), _ddx(grid, ix+1, iq2, id) * dx);
    double vh = _interpolateCubic(tx, grid.xf(ix,   iq2+1, id), _ddx(grid, ix, iq2+1, id) * dx,
				  grid.xf(ix+1, iq2+1, id), _ddx(grid, ix+1, iq2+1, id) * dx);
    
    // Derivatives in Q2
    double vdl, vdh;
    if (q2_lower) {
      // Forward difference for lower q
      vdl = (vh - vl) / dq_1;
      // Central difference for higher q
      double vhh = _interpolateCubic(tx, grid.xf(ix, iq2+2, id), _ddx(grid, ix, iq2+2, id) * dx,
				     grid.xf(ix+1, iq2+2, id), _ddx(grid, ix+1, iq2+2, id) * dx);
      vdh = (vdl + (vhh - vh)/dq_2) / 2.0;
    }
    else if (q2_upper) {
      // Backward difference for higher q
      vdh = (vh - vl) / dq_1;
      // Central difference for lower q
      double vll = _interpolateCubic(tx, grid.xf(ix, iq2-1, id), _ddx(grid, ix, iq2-1, id) * dx,
				     grid.xf(ix+1, iq2-1, id), _ddx(grid, ix+1, iq2-1, id) * dx);
      vdl = (vdh + (vl - vll)/dq_0) / 2.0;
    }
    else {
      // Central difference for both q
      double vll = _interpolateCubic(tx, grid.xf(ix, iq2-1, id), _ddx(grid, ix, iq2-1, id) * dx,
				     grid.xf(ix+1, iq2-1, id), _ddx(grid, ix+1, iq2-1, id) * dx);
      vdl = ( (vh - vl)/dq_1 + (vl - vll)/dq_0 ) / 2.0;
      double vhh = _interpolateCubic(tx, grid.xf(ix, iq2+2, id), _ddx(grid, ix, iq2+2, id) * dx,
				     grid.xf(ix+1, iq2+2, id), _ddx(grid, ix+1, iq2+2, id) * dx);
      vdh = ( (vh - vl)/dq_1 + (vhh - vh)/dq_2 ) / 2.0;
    }

    vdl *= dq;
    vdh *= dq;

    return _interpolateCubic(tq, vl, vdl, vh, vdh);
  }


}
