#pragma once

#include "LHAPDF/types.h"

namespace LHAPDF {


  // Forward declaration
  class PDFGrid;


  /// The general interface for interpolating between grid points.
  class Interpolator {
  public:

    /// Destructor to allow inheritance
    virtual ~Interpolator() {};

    /// Bind to PDFGrid.
    ///
    /// \param pdf the PDFGrid to bind
    virtual void bind(const PDFGrid& pdf) { }

    /// Unbind from PDFGrid.
    ///
    /// \param pdf the PDFGrid to unbind
    virtual void unbind(const PDFGrid& pdf) { }

    /// Interpolate a single-point in (x,Q).
    ///
    /// \param pdf the PDFGrid to interpolate
    /// \param id the MC Parton ID
    /// \param x the momentum fraction
    /// \param q the energy scale
    /// \return
    virtual double interpolateQ(const PDFGrid& pdf, PID_t id, double x, double q) const {
      return interpolateQ2(pdf, id, x, q*q);
    }

    /// Interpolate single-point in Q2.
    ///
    /// \param pdf the PDFGrid to interpolate
    /// \param id the MC Parton ID
    /// \param x the momentum fraction
    /// \param q2 the squared energy scale
    /// \return
    virtual double interpolateQ2(const PDFGrid& pdf, PID_t id, double x, double q2) const = 0;

  };


}
