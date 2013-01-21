#pragma once

#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Forward declaration
  class GridPDF;


  /// The general interface for interpolating between grid points
  class Interpolator {
  public:

    /// Destructor to allow inheritance
    virtual ~Interpolator() { }

    /// Bind to a GridPDF
    void bind(const GridPDF* pdf) { _pdf = pdf; }

    /// Unbind from GridPDF
    void unbind() { _pdf = 0; }

    /// Identify whether this Interpolator has an associated PDF
    bool hasPdf() { return _pdf != 0; }

    /// Get the associated GridPDF
    const GridPDF& pdf() const { return *_pdf; }

    /// Interpolate a single-point in (x,Q)
    virtual double interpolateXQ(int id, double x, double q) const {
      return interpolateXQ2(id, x, q*q);
    }

    /// Interpolate single-point in Q2
    virtual double interpolateXQ2(int id, double x, double q2) const = 0;

    /// @todo The all-PID version of interpolateQ and Q2


  private:

    const GridPDF* _pdf;

  };


}
