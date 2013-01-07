#pragma once

#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Forward declaration
  class PDFGrid;


  /// The general interface for extrapolating beyond grid boundaries.
  class Extrapolator {
  public:

    /// Destructor to allow inheritance
    virtual ~Extrapolator() { }

    /// Bind to a PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    void bind(const PDFGrid* pdf) { _pdf = pdf; }

    /// Unbind from PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    void unbind() { _pdf = 0; }

    /// Get the associated PDFGrid.
    /// @todo Clarify the ownership: who's responsible for deleting? Singleton status?
    const PDFGrid* pdf() const { return _pdf; }

    /// Extrapolate a single-point in (x,Q).
    double extrapolateXQ(int id, double x, double q) const {
      return extrapolateXQ2(id, x, q*q );
    }

    /// Extrapolate a single-point in (x,Q2).
    virtual double extrapolateXQ2(int, double x, double q2) const = 0;

    /// @todo The all-PID version of extrapolateQ and Q2


  private:

    const PDFGrid* _pdf;

  };


}
