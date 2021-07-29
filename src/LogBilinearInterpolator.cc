// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/LogBilinearInterpolator.h"

namespace LHAPDF {


  namespace { // Unnamed namespace

    // One-dimensional linear interpolation for y(x)
    inline double _interpolateLinear(double x, double xl, double xh, double yl, double yh)	{
      assert(x >= xl);
      assert(xh >= x);
      return yl + (x - xl) / (xh - xl) * (yh - yl);
    }

  }

  
  double LogBilinearInterpolator::_interpolateXQ2(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, int id) const {
    if (grid.xsize() < 2)
      throw GridError("PDF subgrids are required to have at least 2 x-knots for use with LogBilinearInterpolator");
    if (grid.q2size() < 2)
      throw GridError("PDF subgrids are required to have at least 2 Q2-knots for use with LogBilinearInterpolator");

    // manual cache some values
    const double logx = log(x);
    const double logx0 = grid.logxs(ix);
    const double logx1 = grid.logxs(ix+1);
    
    // First interpolate in x
    const double f_ql = _interpolateLinear(logx, logx0, logx1, grid.xf(ix, iq2, id), grid.xf(ix+1, iq2, id));
    const double f_qh = _interpolateLinear(logx, logx0, logx1, grid.xf(ix, iq2+1, id), grid.xf(ix+1, iq2+1, id));
    // Then interpolate in Q2, using the x-ipol results as anchor points
    return _interpolateLinear(log(q2), grid.logq2s(iq2), grid.logq2s(iq2+1), f_ql, f_qh);
  }


}
