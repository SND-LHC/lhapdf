// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Extrapolator.h"

namespace LHAPDF {


  /// Extrapolates using the closest point on the Grid.
  class NearestPointExtrapolator : public Extrapolator {
  public:

    double extrapolateXQ2(int id, double x, double q2) const;

  };


}
