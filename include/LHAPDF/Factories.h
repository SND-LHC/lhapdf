#pragma once

#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/BilinearInterpolator.h"
#include "LHAPDF/BicubicInterpolator.h"
#include "LHAPDF/ErrExtrapolator.h"
#include "LHAPDF/NearestPointExtrapolator.h"

namespace LHAPDF {


  /// Create a new PDF from the given data file path.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  inline PDF* mkPDF(const std::string& path) {
    // First create an Info object to work out what format of PDF this is:
    LHAPDF::Info info(path);
    const string fmt = info.metadata("Format");
    // Then use the format information to call the appropriate concrete PDF constructor:
    if (fmt == "lhagrid") return new PDFGrid(path);
    throw FactoryError("No LHAPDF factory defined for format type '" + fmt + "'");
  }

  /// Create a new PDF with the given PDF set name and member ID.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  inline PDF* mkPDF(const std::string& setname, int member) {
    path mempath = findFile(pdfmempath(setname, member));
    return mkPDF(mempath.native());
  }

  /// Create a new PDF with the given LHAPDF ID code.
  ///
  /// Returns a 'new'ed PDF by pointer.
  /// The caller is responsible for deletion of the created object.
  inline PDF* mkPDF(int lhaid) {
    const pair<string,int> setname_memid = lookupPDF(lhaid);
    return mkPDF(setname_memid.first, setname_memid.second);
  }


  /// Interpolator factory
  ///
  /// Returns a 'new'ed Interpolator by pointer. Unless passed to a PDFGrid,
  /// the caller is responsible for deletion of the created object.
  inline Interpolator* mkInterpolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = to_lower_copy(name);
    if (iname == "linear")
      return new BilinearInterpolator();
    else if (iname == "cubic")
      return new BicubicInterpolator();
    else
      throw FactoryError("Undeclared interpolator requested: " + name);
  }


  /// Extrapolator factory
  ///
  /// Returns a 'new'ed Extrapolator by pointer. Unless passed to a PDFGrid,
  /// the caller is responsible for deletion of the created object.
  inline Extrapolator* mkExtrapolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = to_lower_copy(name);
    if (iname == "nearest")
      return new NearestPointExtrapolator();
    else if (iname == "error")
      return new ErrExtrapolator();
    else
      throw FactoryError("Undeclared extrapolator requested: " + name);
  }


  // /// AlphaS factory
  // ///
  // /// Returns a 'new'ed AlphaS by pointer. Unless attached to a PDF,
  // /// the caller is responsible for deletion of the created object.
  // inline AlphaS* mkAlphaS(const std::string& name) {
  //   // Convert name to lower case for comparisons
  //   const std::string iname = to_lower_copy(name);
  //   if (iname == "analytic") //< Need a better name (and class name) if we implement more than one analytic approximation
  //     return new AlphaS_Analytic();
  //   else if (iname == "ode")
  //     return new AlphaS_ODE();
  //   else if (iname == "ipol")
  //     return new AlphaS_Ipol();
  //   else
  //     throw FactoryError("Undeclared AlphaS requested: " + name);
  // }


}
