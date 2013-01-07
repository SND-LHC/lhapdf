#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/ErrExtrapolator.h"
#include <iostream>

namespace LHAPDF {


  double ErrExtrapolator::extrapolateQ2(const PDFGrid& grid, int id, double x, double q2) const {
    std::cerr << "Extrapolation not allowed." << std::endl;
    throw -1;
  }


}
