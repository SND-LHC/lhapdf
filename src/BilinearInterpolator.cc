#include "LHAPDF/GridPDF.h"
#include "LHAPDF/BilinearInterpolator.h"
#include <algorithm>
#include <cassert>

using namespace std;

namespace LHAPDF {


  namespace { // Unnamed namespace

    // One-dimensional linear interpolation for y(x)
    inline double _interpolateLinear(double x, double xl, double xh, double yl, double yh)	{
      assert(x >= xl);
      assert(xh >= x);
      return yl + (x - xl) / (xh - xl) * (yh - yl);
    }

  }


  double BilinearInterpolator::interpolateXQ2(int id, double x, double q2) const {
    /// @todo Move the following to the Interpolator interface and implement caching
    // Subgrid and index lookup
    const GridPDF::KnotArray1F& subgrid = pdf().subgrid(id, q2);
    const size_t ix = subgrid.xlow(x);
    const size_t iq2 = subgrid.q2low(q2);

    // First interpolate in x
    const double f_ql = _interpolateLinear(x, subgrid.xs()[ix], subgrid.xs()[ix+1], subgrid.xf(ix, iq2), subgrid.xf(ix+1, iq2));
    const double f_qh = _interpolateLinear(x, subgrid.xs()[ix], subgrid.xs()[ix+1], subgrid.xf(ix, iq2+1), subgrid.xf(ix+1, iq2+1));

    // Then interpolate in Q2, using the x-ipol results as anchor points
    return _interpolateLinear(q2, subgrid.q2s()[iq2], subgrid.q2s()[iq2+1], f_ql, f_qh);
  }


}
