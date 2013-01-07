#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// @brief Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class BicubicInterpolator : public Interpolator {
  public:
    double interpolateXQ2(int id, double x, double q2) const;

  private:
    double ddx(const double* xfs, size_t xidx, size_t q2idx) const;
  };


  /// @todo Make a log-space BicubicLogInterpolator based on this

}
