#pragma once

#include "LHAPDF/Extrapolator.h"

namespace LHAPDF {


  /// Extrapolates using the closest point on the Grid.
  class NearestPointExtrapolator : public Extrapolator {

    /// \param id The PDG parton ID code
    /// \param x The momentum fraction
    /// \param q2 The squared energy scale
    /// \return The xf value at (x,q2)
    double extrapolateXQ2(int id, double x, double q2) const;

  };


}
