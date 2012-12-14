#pragma once

#include "LHAPDF/types.h"

namespace LHAPDF {


  // Forward declaration
  class PDFGrid;


  /// The general interface for extrapolating beyond grid boundaries.
  class Extrapolator {
  public:

    /// Destructor to allow inheritance
    virtual ~Extrapolator() { }


    /// Bind to PDFGrid.
    ///
    /// \param pdf the PDFGrid to bind
    virtual void bind(const PDFGrid& pdf) { }


    /// Unbind from PDFGrid.
    ///
    /// \param pdf the PDFGrid to unbind
    virtual void unbind( const PDFGrid& pdf) { }


    /// Extrapolate a single-point in (x,Q).
    ///
    /// \param Grid the PDFGrid to extrapolate
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q the energy scale
    /// \return
    double extrapolateQ(const PDFGrid& pdf, PID_t id, double x, double q) const {
      return extrapolateQ2( pdf, id, x, q*q );
    }

    /// Extrapolate a single-point in (x,Q2).
    ///
    /// \param Grid the PDFGrid to extrapolate
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    virtual double extrapolateQ2( const PDFGrid&, PID_t, double x, double q2 ) const = 0;

    /// @todo The all-PID version of extrapolateQ and Q2

  };


}
