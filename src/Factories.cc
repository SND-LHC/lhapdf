// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/BilinearInterpolator.h"
#include "LHAPDF/BicubicInterpolator.h"
#include "LHAPDF/LogBilinearInterpolator.h"
#include "LHAPDF/LogBicubicInterpolator.h"
#include "LHAPDF/ErrExtrapolator.h"
#include "LHAPDF/NearestPointExtrapolator.h"
#include "LHAPDF/AlphaS.h"

namespace LHAPDF {


  Info& getConfig() {
    return Config::get();
  }


  PDFSet& getPDFSet(const string& setname) {
    static map<string, PDFSet> _sets;
    map<string, PDFSet>::iterator it = _sets.find(setname);
    if (it != _sets.end()) return it->second;
    _sets[setname] = PDFSet(setname);
    return _sets[setname];
  }


  PDF* mkPDF(const string& setname, int member) {
    // First create an Info object to work out what format of PDF this is:
    Info info(findpdfmempath(setname, member));
    const string fmt = info.metadata("Format");
    // Then use the format information to call the appropriate concrete PDF constructor:
    if (fmt == "lhagrid1") return new GridPDF(setname, member);
    /// @todo Throw a deprecation error if format version is too old or new
    throw FactoryError("No LHAPDF factory defined for format type '" + fmt + "'");
  }

  PDF* mkPDF(int lhaid) {
    const pair<string,int> setname_memid = lookupPDF(lhaid);
    return mkPDF(setname_memid.first, setname_memid.second);
  }


  void mkPDFs(const string& setname, vector<PDF*>& pdfs) {
    getPDFSet(setname).mkPDFs(pdfs);
  }

  vector<PDF*> mkPDFs(const string& setname) {
    return getPDFSet(setname).mkPDFs();
  }


  Interpolator* mkInterpolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = to_lower_copy(name);
    if (iname == "linear")
      return new BilinearInterpolator();
    else if (iname == "cubic")
      return new BicubicInterpolator();
    else if (iname == "log")
      return new LogBilinearInterpolator();
    else if (iname == "logcubic")
      return new LogBicubicInterpolator();
    else
      throw FactoryError("Undeclared interpolator requested: " + name);
  }


  Extrapolator* mkExtrapolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = to_lower_copy(name);
    if (iname == "nearest")
      return new NearestPointExtrapolator();
    else if (iname == "error")
      return new ErrExtrapolator();
    else
      throw FactoryError("Undeclared extrapolator requested: " + name);
  }


  AlphaS* mkAlphaS(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = to_lower_copy(name);
    if (iname == "analytic") //< Need a better name (and class name) if we implement more than one analytic approximation
      return new AlphaS_Analytic();
    else if (iname == "ode")
      return new AlphaS_ODE();
    /// @todo Reinstate
    // else if (iname == "ipol")
    //   return new AlphaS_Ipol();
    else
      throw FactoryError("Undeclared AlphaS requested: " + name);
  }


}
