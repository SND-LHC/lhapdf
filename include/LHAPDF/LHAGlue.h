// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

/// @file LHAGlue.h
/// A file that provides backwards compatibility for some C functions from lhapdf 5.x

namespace LHAPDF {
  /// Number of members available in the current set.
  int numberPDF();
  /// Number of members available in the current set.
  int numberPDF(int nset);
  /// The choice of PDF member out of one distribution.
  void initPDF(int memset);
  /// The choice of PDF member out of one distribution.
  void initPDF(int nset, int memset);
  /// Nucleon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  double xfx(double x, double Q, int fl);
  /// Nucleon PDF: returns @c x f(x, Q) for flavour @a fl - flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  double xfx(int nset, double x, double Q, int fl);


  /// @brief Distinction between evolution or interpolation PDF sets.
  /// Enum to choose whether evolution (i.e. @c LHpdf data file) or
  /// interpolation (i.e. @c LHgrid data file) is used.
  enum SetType {
    EVOLVE = 0, LHPDF = 0,
    INTERPOLATE = 1, LHGRID = 1
  };



}

