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
    /// @todo Move to be "normal" functions
    double interpolateXQ2(const GridPDF::KnotArray1F& subgrid, size_t ix, size_t iq2) const;
    double ddx(const GridPDF::KnotArray1F& subgrid, size_t ix, size_t iq2) const;
  };

  /// @todo Make a log-space BicubicLogInterpolator based on this


}
