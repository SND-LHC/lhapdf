#pragma once

namespace LHAPDF {


  // Forward declarations to avoid circular dependencies
  class PDF;
  class Interpolator;
  class Extrapolator;
  class AlphaS;


  /// @name Factory functions for making PDFs and other objects
  //@{

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
  /// Returns a 'new'ed Interpolator by pointer. Unless passed to a GridPDF,
  /// the caller is responsible for deletion of the created object.
  Interpolator* mkInterpolator(const std::string& name);


  /// Extrapolator factory
  ///
  /// Returns a 'new'ed Extrapolator by pointer. Unless passed to a GridPDF,
  /// the caller is responsible for deletion of the created object.
  Extrapolator* mkExtrapolator(const std::string& name);


  /// AlphaS factory
  ///
  /// Returns a 'new'ed AlphaS by pointer. Unless attached to a PDF,
  /// the caller is responsible for deletion of the created object.
  AlphaS* mkAlphaS(const std::string& name);

  //@}


}
