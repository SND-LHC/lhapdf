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

    /// @brief x-caching struct definition
    struct XCache {
      // Subgrid hash
      size_t xhash = 0;
      // Defining params from last call (initialised to unphysical values, so first use will set the cache)
      double x = -1;
      size_t ix = SIZE_MAX;
      // Cached params
      double logx;
      double dlogx_1;
      double tlogx;
    };

    /// @brief Q2-caching struct definition
    struct Q2Cache {
      // Subgrid hashes
      size_t q2hash = 0;
      // Defining params from last call (initialised to unphysical values, so first use will set the cache)
      double q2 = -1;
      size_t iq2 = SIZE_MAX;
      // Cached params
      double logq2;
      double dlogq_0;
      double dlogq_1;
      double dlogq_2;
      double tlogq;
    };

    /// @brief Get and update the current caching structs for interpolation params
    ///
    /// @note Caches are handled separately for x and Q since they can be sampled very differently.
    ///@{
    using XCacheMap = map<size_t,XCache>;
    static XCache& _getCacheX(const KnotArray1F& subgrid, double x, size_t ix);

    using Q2CacheMap = map<size_t,Q2Cache>;
    static Q2Cache& _getCacheQ2(const KnotArray1F& subgrid, double q2, size_t iq2);
    ///@}

  };


}

#endif
