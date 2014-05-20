// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/PDF.h"

namespace LHAPDF {


  /// @name Single beam reweighting functions
  //@{

  /// Get the PDF reweighting factor for a beam with id,x,Q parameters, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  inline double weightxQ2(int id, double x, double Q2, const PDF& basepdf, const PDF& newpdf) {
    const double xfx_base = basepdf.xfxQ2(id, x, Q2);
    const double xfx_new = newpdf.xfxQ2(id, x, Q2);
    return xfx_new / xfx_base;
  }

  /// Get the PDF reweighting factor for a beam with id,x,Q parameters, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  template <typename PDFPTR>
  inline double weightxQ2(int id, double x, double Q2, const PDFPTR basepdf, const PDFPTR newpdf) {
    return weightxQ2(id, x, Q2, *basepdf, *newpdf);
  }

  /// Get the PDF reweighting factor for a beam with id,x,Q parameters, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  inline double weightxQ(int id, double x, double Q, const PDF& basepdf, const PDF& newpdf) {
    return weightxQ2(id, x, sqr(Q), basepdf, newpdf);
  }

  /// Get the PDF reweighting factor for a beam with id,x,Q parameters, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  template <typename PDFPTR>
  inline double weightxQ(int id, double x, double Q, const PDFPTR basepdf, const PDFPTR newpdf) {
    return weightxQ(id, x, Q, *basepdf, *newpdf);
  }

  //@}


  /// @name Two-beam reweighting functions
  //@{

  /// Get the PDF reweighting factor for two beams, one with id1,x1 and the other with id2,x2, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  inline double weightxxQ2(int id1, int id2, double x1, double x2, double Q2, const PDF& basepdf, const PDF& newpdf) {
    const double w1 = weightxQ2(id1, x1, Q2, basepdf, newpdf);
    const double w2 = weightxQ2(id2, x2, Q2, basepdf, newpdf);
    return w1 * w2;
  }

  /// Get the PDF reweighting factor for two beams, one with id1,x1 and the other with id2,x2, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  template <typename PDFPTR>
  inline double weightxxQ2(int id1, int id2, double x1, double x2, double Q2, const PDFPTR basepdf, const PDFPTR newpdf) {
    return weightxxQ2(id1, id2, x1, x2, Q2, *basepdf, *newpdf);
  }

  /// Get the PDF reweighting factor for two beams, one with id1,x1 and the other with id2,x2, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  inline double weightxxQ(int id1, int id2, double x1, double x2, double Q, const PDF& basepdf, const PDF& newpdf) {
    return weightxxQ2(id1, id2, x1, x2, sqr(Q), basepdf, newpdf);
  }

  /// Get the PDF reweighting factor for two beams, one with id1,x1 and the other with id2,x2, from basepdf to newpdf
  /// @note For NLO calculations, in general different PDF values enter for each counterterm: be careful.
  template <typename PDFPTR>
  inline double weightxxQ(int id1, int id2, double x1, double x2, double Q, const PDFPTR basepdf, const PDFPTR newpdf) {
    return weightxxQ(id1, id2, x1, x2, Q, *basepdf, *newpdf);
  }

  //@}


}
