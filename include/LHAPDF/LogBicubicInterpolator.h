// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_LogBicubicInterpolator_H
#define LHAPDF_LogBicubicInterpolator_H

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// @brief Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class LogBicubicInterpolator : public Interpolator {
  public:

    /// Implementation of (x,Q2) interpolation
    double _interpolateXQ2(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) const;

    /// @brief Caching struct definition
    ///
    /// @todo Cache more ipol-weight variables?
    struct XQ2Cache {
      // Subgrid hashes
      size_t xhash = 0;
      size_t q2hash = 0;
      // Defining params from last call (initialised to unphysical values, so first use will set the cache)
      double x = -1;
      size_t ix = 0;
      double q2 = -1;
      size_t iq2 = 0;
      // Cached params
      double logx;
      double logq2;
      double dlogx_1;
      double tlogx;
      double dlogq_0;
      double dlogq_1;
      double dlogq_2;
      double tlogq;
    };

    /// @brief Get and update the current caching struct for interpolation params
    ///
    /// @note Implemented as a thread-safe Meyers Singleton.
    /// @note Cache is handled separately for x and Q since they can be sampled very differently.
    static XQ2Cache& _getCache(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2);

  };


}

#endif
