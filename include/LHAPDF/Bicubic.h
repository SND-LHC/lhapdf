#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class BicubicInterpolator : public Interpolator {
  public:

    /// \param Grid the PDFGrid
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    double interpolateQ2(const PDFGrid& pdf, PID_t id, double x, double q2) const;

  };


}
