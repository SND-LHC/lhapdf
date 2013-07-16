// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

/// @file LHAGlue.h
/// A file that provides backwards compatibility for some C functions from LHAPDF 5.x


/// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
/// @todo Why defined outside the LHAPDF namespace?!
/// @deprecated Use the proper C++ interface of LHAPDF6 instead!
void initPDFSetByName(const std::string& filename);

/// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
/// @todo Why defined outside the LHAPDF namespace?!
/// @deprecated Use the proper C++ interface of LHAPDF6 instead!
void initPDFSetByName(int nset, const std::string& filename);

namespace LHAPDF {

  /// Number of members available in the current set.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int numberPDF();

  /// Number of members available in the current set.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int numberPDF(int nset);

  /// The choice of PDF member out of one distribution.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDF(int memset);

  /// The choice of PDF member out of one distribution.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDF(int nset, int memset);

  /// Nucleon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - flavour encoding as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double xfx(double x, double Q, int fl);

  /// Nucleon PDF: returns @c x f(x, Q) for flavour @a fl - flavour encoding as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double xfx(int nset, double x, double Q, int fl);

  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void xfx(double x, double Q, double* results);

  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void xfx(int nset, double x, double Q, double* results);

  /// Nucleon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  std::vector<double> xfx(double x, double Q);

  /// Nucleon PDF: returns a vector @c x f_i(x, Q) with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  std::vector<double> xfx(int nset, double x, double Q);

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double alphasPDF(double Q);

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double alphasPDF(int nset, double Q);

  /// @brief Use @a member in current PDF set.
  /// This operation is computationally cheap.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void usePDFMember(int member);

  /// @brief Use @a member in PDF set @a nset (multi-set version).
  /// This operation is computationally cheap.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void usePDFMember(int nset, int member);

  /// Initialise @a member in PDF set @a setid.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDFSet(const std::string& name);

}
