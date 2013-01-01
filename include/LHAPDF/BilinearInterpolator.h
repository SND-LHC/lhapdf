#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// Implementation of bilinear interpolation
  class BilinearInterpolator : public Interpolator {
  public:
    double interpolateXQ2(PID_t id, double x,  double q2) const;
  };


  /// @todo Make a log-space BilinearLogInterpolator based on this


}
