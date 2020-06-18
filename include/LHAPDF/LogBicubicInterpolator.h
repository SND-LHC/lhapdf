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

    /// @brief A single set of cached x-variables
    struct XCache {
      /// Defining params from call (initialised to unphysical values, so first use will set the cache)
      double x = -1;
      //size_t ix;

      /// Cached params
      double logx;
      double dlogx_1;
      double tlogx;
    };

    /// A multi-level x-variable cache for a single subgrid hash
    struct XCaches {
      // Number of cache levels
      static const size_t N = 4;

      // Latest-call index
      size_t ilast = 0;

      // List of N cached-value sets
      array<XCache, N> caches;

      XCache& operator[] (size_t index) { return caches[index]; }

    };


    /// @brief A single set of cached Q2-variables
    struct Q2Cache {
      /// Defining params from call (initialised to unphysical values, so first use will set the cache)
      double q2 = -1;

      /// Cached params
      double logq2;
      double dlogq_0;
      double dlogq_1;
      double dlogq_2;
      double tlogq;
    };

    /// A multi-level Q2-variable cache for a single subgrid hash
    struct Q2Caches {
      // Number of cache levels
      static const size_t N = 4;

      // Latest-call index
      size_t ilast = 0;

      // List of N cached-value sets
      array<Q2Cache, N> caches;

      Q2Cache& operator[] (size_t index) { return caches[index]; }

    };


    /// @brief Get and update the current caching structs for interpolation params
    ///
    /// @note Caches are handled separately for x and Q since they can be sampled very differently.
    ///@{
    using XCachesMap = map<size_t,XCaches>;
    static XCache& _getCacheX(const KnotArray1F& subgrid, double x, size_t ix);

    using Q2CachesMap = map<size_t,Q2Caches>;
    static Q2Cache& _getCacheQ2(const KnotArray1F& subgrid, double q2, size_t iq2);
    ///@}

  };


}

#endif
