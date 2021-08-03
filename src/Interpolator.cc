// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/GridPDF.h"

namespace LHAPDF {


  double Interpolator::interpolateXQ2(int id, double x, double q2) const {
    // Subgrid lookup
    /// @todo Do this in two stages to cache the KnotArrayNF?
    /// @todo Add flavour error checking

    const KnotArray& grid = pdf().knotarray();
    /// @todo Cache this index lookup for performance?
    //  Maybe compiler finds this, as this is now *always* the same call

    // 25% of computing time spend here
    const size_t ix  = grid.ixbelow(x);
    const size_t iq2 = grid.iq2below(q2);
    
    /// Call the overloaded interpolation routine on this subgrid
    return _interpolateXQ2(grid, x, ix, q2, iq2, id);
  }

  void Interpolator::interpolateXQ2(double x, double q2, std::vector<double>& ret) const {
    // Subgrid lookup
    /// @todo Do this in two stages to cache the KnotArrayNF?
    /// @todo Add flavour error checking

    const KnotArray& grid = pdf().knotarray();
    /// @todo Cache this index lookup for performance?
    //  Maybe compiler finds this, as this is now *always* the same call

    // 25% of computing time spend here
    const size_t ix  = grid.ixbelow(x);
    const size_t iq2 = grid.iq2below(q2);
    
    /// Call the overloaded interpolation routine on this subgrid
    _interpolateXQ2(grid, x, ix, q2, iq2, ret);
  }
  

}
