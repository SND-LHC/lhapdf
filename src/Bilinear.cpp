#include "Bilinear.h"
#include "PDFGrid.h"

/// The one dimensional implementation of linear interpolation.

/// \param X
/// \param XL
/// \param XH
/// \param YL
/// \param YH
/// \return 
inline double interpolateLinear( const double X, const double XL, const double XH, const double YL, const double YH )
{
    return YL + ( X - XL ) / ( XH - XL ) * ( YH - YL );
}

double BilinearInterpolator::interpolateQ2( const PDFGrid& PDF, const PID_t ID, const X_t X, const Q2_t Q2 ) const
{
    //Find top left corner of patch
    Idx_t idx_x, idx_q2;
    PDF.lookupClosestLow( X, Q2, idx_x, idx_q2 );
    
    //Calculate interpolation along X axis
    const AxisKnots xpoints = PDF.getXKnots();
    const AxisKnots qpoints = PDF.getQ2Knots();
    const PIDdata pid = PDF.getPIDData( ID );
    
    double f_ql = interpolateLinear( X, xpoints[idx_x], xpoints[idx_x+1], pid[PDF.index(idx_x, idx_q2)], pid[PDF.index(idx_x+1, idx_q2)] );
    double f_qh = interpolateLinear( X, xpoints[idx_x], xpoints[idx_x+1], pid[PDF.index(idx_x, idx_q2+1)], pid[PDF.index(idx_x+1, idx_q2+1)] );
    
    //Calculate interpolation along Q axis, using both values of X axis
    return interpolateLinear( Q2, qpoints[idx_q2], qpoints[idx_q2+1], f_ql, f_qh );
}
