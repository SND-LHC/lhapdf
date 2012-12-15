#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/BilinearInterpolator.h"

namespace LHAPDF {


  /// Two-dimensional linear interpolation.
  ///
  /// \param X
  /// \param XL
  /// \param XH
  /// \param YL
  /// \param YH
  /// \return
  inline double interpolateLinear(double X, double XL, double XH, double YL, double YH )	{
    return YL + ( X - XL ) / ( XH - XL ) * ( YH - YL );
  }


  double BilinearInterpolator::interpolateQ2(const PDFGrid& grid, PID_t id, double x, double q2) const {
    // Find top left corner of patch
    size_t idx_x, idx_q2;
    grid.lookupClosestLow( x, q2, idx_x, idx_q2 );

    // Calculate interpolation along X axis
    const AxisKnots xpoints = grid.getXKnots();
    const AxisKnots qpoints = grid.getQ2Knots();
    const PIDdata pid = grid.getPIDData( id );

    double f_ql = interpolateLinear( x, xpoints[idx_x], xpoints[idx_x+1], pid[grid.index(idx_x, idx_q2)], pid[grid.index(idx_x+1, idx_q2)] );
    double f_qh = interpolateLinear( x, xpoints[idx_x], xpoints[idx_x+1], pid[grid.index(idx_x, idx_q2+1)], pid[grid.index(idx_x+1, idx_q2+1)] );

    // Calculate interpolation along Q axis, using both values of X axis
    return interpolateLinear( q2, qpoints[idx_q2], qpoints[idx_q2+1], f_ql, f_qh );
  }


}
