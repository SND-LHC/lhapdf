#include "LHAPDF/Nearest.h"
#include "LHAPDF/PDFGrid.h"

namespace LHAPDF {


  double NearestExtrapolator::extrapolateQ2( const PDFGrid& grid, PID_t id, const X_t x, const Q2_t q2 ) const {
    Idx_t xidx, q2idx;
    grid.lookupClosest( x, q2, xidx, q2idx );
    return grid.getPIDData(id)[grid.index(xidx, q2idx)];
  }


}
