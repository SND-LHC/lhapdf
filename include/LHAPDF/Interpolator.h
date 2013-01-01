#pragma once

#include "LHAPDF/Types.h"

namespace LHAPDF {


  // Forward declaration
  class PDFGrid;


  /// The general interface for interpolating between grid points.
  class Interpolator {
  public:

    /// Destructor to allow inheritance
    virtual ~Interpolator() { }

    /// Bind to a PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    void bind(const PDFGrid* pdf) { _pdf = pdf; }

    /// Unbind from PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    void unbind() { _pdf = 0; }

    /// Get the associated PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    const PDFGrid* pdf() const { return _pdf; }

    /// Interpolate a single-point in (x,Q).
    virtual double interpolateXQ(PID_t id, double x, double q) const {
      return interpolateXQ2(id, x, q*q);
    }

    /// Interpolate single-point in Q2.
    virtual double interpolateXQ2(PID_t id, double x, double q2) const = 0;


  private:

    const PDFGrid* _pdf;

  };


}
