#include "LHAPDF/NearestPointExtrapolator.h"
#include "LHAPDF/PDFGrid.h"

namespace LHAPDF {


  double NearestPointExtrapolator::extrapolateXQ2(int id, double x, double q2) const {
    /// @todo Functionality needs to be re-added to PDFGrid... or implemented here.
    // pair<size_t, size_t> idx = pdf()->lookupClosest(x, q2);
    // xidx = idx.first;
    // q2idx = idx.second;
    // return pdf()->getPIDData(id)[pdf()->ptindex(xidx, q2idx)];
    return 0;
  }


}
