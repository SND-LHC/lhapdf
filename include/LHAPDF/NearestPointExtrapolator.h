#pragma once

#include "LHAPDF/Extrapolator.h"

namespace LHAPDF {


  /// Extrapolates using the closest point on the Grid.
  class NearestPointExtrapolator : public Extrapolator {

    /// \param Grid the GridPDF to extrapolate
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return the xfx value at (X,Q2)
    double extrapolateXQ2(int, double x, double q2 ) const;

  };


}
