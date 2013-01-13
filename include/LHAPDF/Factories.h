#pragma once

#include "LHAPDF/PDF.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/Extrapolator.h"
#include "LHAPDF/AlphaS.h"

namespace LHAPDF {


  /// Create a new PDF from the given data file path.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  PDF* mkPDF(const std::string& path);

  /// Create a new PDF with the given PDF set name and member ID.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  PDF* mkPDF(const std::string& setname, int member);

  /// Create a new PDF with the given LHAPDF ID code.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  PDF* mkPDF(int lhaid);


  /// Interpolator factory
  ///
  /// Returns a 'new'ed Interpolator by pointer. Unless passed to a PDFGrid,
  /// the caller is responsible for deletion of the created object.
  Interpolator* mkInterpolator(const std::string& name);


  /// Extrapolator factory
  ///
  /// Returns a 'new'ed Extrapolator by pointer. Unless passed to a PDFGrid,
  /// the caller is responsible for deletion of the created object.
  Extrapolator* mkExtrapolator(const std::string& name);


  /// AlphaS factory
  ///
  /// Returns a 'new'ed AlphaS by pointer. Unless attached to a PDF,
  /// the caller is responsible for deletion of the created object.
  AlphaS* mkAlphaS(const std::string& name);


}
