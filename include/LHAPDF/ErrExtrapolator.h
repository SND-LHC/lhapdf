#pragma once

#include "LHAPDF/Extrapolator.h"

namespace LHAPDF {

  /// Extrapolates using the closest point on the Grid.
  class ErrExtrapolator : public Extrapolator {
  public:

    double extrapolateXQ2(int id, double x, double q2) const;

  };


}
