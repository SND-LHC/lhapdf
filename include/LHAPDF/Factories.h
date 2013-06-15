// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include <string>

namespace LHAPDF {


  // Forward declarations to avoid circular dependencies
  class PDF;
  class PDFSet;
  class Config;
  class Interpolator;
  class Extrapolator;
  class AlphaS;


  /// Get the global configuration object
  ///
  /// The global config is populated by reading from lhapdf.conf if it is found
  /// in the search paths. It is a singleton, hence the 'get' rather than 'mk'
  /// function name.
  ///
  /// @note The LHAPDF system is responsible for deletion of the returned
  /// object. Do NOT delete it yourself! Hence the return by reference rather
  /// than pointer.
  // Config& getConfig();
  Info& getConfig();


  /// Get the PDFSet with the given set name.
  ///
  /// Returns a PDFSet by reference. When this function is used for
  /// access, only one PDFSet object is made per set name... hence the
  /// 'get' rather than 'mk' function name.
  ///
  /// This function is intended particularly for use where it would be
  /// inefficient to have to repeatedly construct a PDFSet by name. The
  /// canonical use case is internal: the Info system uses this to ensure that
  /// cascading of config settings is efficient, and also allows the automatic
  /// application of set-level changes to all PDF member objects in that set.
  ///
  /// @note The LHAPDF system is responsible for deletion of the returned
  /// object. Do NOT delete it yourself! Hence the return by reference rather
  /// than pointer.
  PDFSet& getPDFSet(const std::string& setname);


  /// @name Factory functions for making single PDF members
  //@{

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

  //@}


  /// @name Factory functions for making all PDF members in a set
  //@{

  /// Get all PDFs in a named set (return by filling the supplied vector).
  void mkPDFs(const std::string& setname, std::vector<PDF*>& pdfs);

  /// Get all PDFs in a named set (return by a new vector).
  std::vector<PDF*> mkPDFs(const std::string& setname);

  /// @todo Add a templated version which can return a vector of smart ptrs

  //@}


  /// @name Factory functions for making other LHAPDF framework objects
  //@{

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
