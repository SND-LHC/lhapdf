#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// @brief Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class BicubicInterpolator : public Interpolator {
  public:
    double interpolateXQ2(const PDFGrid& pdf, PID_t id, double x, double q2) const;
  };


}
