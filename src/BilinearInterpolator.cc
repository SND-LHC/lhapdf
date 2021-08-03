// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/BilinearInterpolator.h"

namespace LHAPDF {


  namespace { // Unnamed namespace

    // One-dimensional linear interpolation for y(x)
    inline double _interpolateLinear(const double x, const double xl, const double xh, const double yl, const double yh){
      assert(x >= xl);
      assert(xh >= x);
      return yl + (x - xl) / (xh - xl) * (yh - yl);
    }

    // Bilinear interpolation for single pid
    double _interpolateSinglePid(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, int id){
      // First interpolate in x
      const double f_ql = _interpolateLinear(x, grid.xs(ix), grid.xs(ix+1), grid.xf(ix, iq2, id), grid.xf(ix+1, iq2, id));
      const double f_qh = _interpolateLinear(x, grid.xs(ix), grid.xs(ix+1), grid.xf(ix, iq2+1, id), grid.xf(ix+1, iq2+1, id));
      // Then interpolate in Q2, using the x-ipol results as anchor points
      return _interpolateLinear(q2, grid.q2s(iq2), grid.q2s(iq2+1), f_ql, f_qh);
    }
  }
  

  double BilinearInterpolator::_interpolateXQ2(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, int id) const {
    if (grid.xsize() < 2)
      throw GridError("PDF subgrids are required to have at least 2 x-knots for use with BilinearInterpolator");
    if (grid.q2size() < 2)
      throw GridError("PDF subgrids are required to have at least 2 Q2-knots for use with BilinearInterpolator");

    return _interpolateSinglePid(grid, x, ix, q2, iq2, id);
  }


  void BilinearInterpolator::_interpolateXQ2(const KnotArray& grid, double x, size_t ix, double q2, size_t iq2, std::vector<double>& ret) const{
    if (grid.xsize() < 2)
      throw GridError("PDF subgrids are required to have at least 2 x-knots for use with BilinearInterpolator");
    if (grid.q2size() < 2)
      throw GridError("PDF subgrids are required to have at least 2 Q2-knots for use with BilinearInterpolator");
    ret.clear();    
    for(int pid(-6); pid <= 6; ++pid){
      int id = grid._lookup[pid + 6];
      if(id == -1){
	ret.push_back(0);
      } else {
	ret.push_back(_interpolateSinglePid(grid, x, ix, q2, iq2, id));
      }
    }
  }

} // End of namespace LHAPDF
