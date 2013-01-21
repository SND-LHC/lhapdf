#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// @brief Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class BicubicInterpolator : public Interpolator {
  public:

    double interpolateXQ2(int id, double x, double q2) const;

  };

  /// @todo Make a log-space BicubicLogInterpolator based on this


}
