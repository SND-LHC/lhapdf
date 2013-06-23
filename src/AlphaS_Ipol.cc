// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {

  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS_Ipol::nf_Q2(double q2) const {
    int nf = 0;
    for (int it = 0; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return nf;
  }

  // Interpolate alpha_s from tabulated points in Q2 via metadata
  double AlphaS_Ipol::alphasQ2(double q2) const {
    /// @todo Implement interpolation-based alpha_s
    return -1;
  }


}
