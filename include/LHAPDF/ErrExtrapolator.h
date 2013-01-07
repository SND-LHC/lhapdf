#pragma once

#include "LHAPDF/Extrapolator.h"

namespace LHAPDF {

  /// Extrapolates using the closest point on the Grid.
  class ErrExtrapolator : public Extrapolator {

    /// \param Grid the PDFGrid to extrapolate
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return the xfx value at (X,Q2)
    double extrapolateQ2(const PDFGrid&, int, double x, double q2) const;

  };


}
