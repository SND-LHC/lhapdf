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


  Interpolator* mkInterpolator(const string& name) {
    // Convert name to lower case for comparisons
    const string iname = to_lower_copy(name);
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


  Extrapolator* mkExtrapolator(const string& name) {
    // Convert name to lower case for comparisons
    const string iname = to_lower_copy(name);
    if (iname == "nearest")
      return new NearestPointExtrapolator();
    else if (iname == "error")
      return new ErrExtrapolator();
    else
      throw FactoryError("Undeclared extrapolator requested: " + name);
  }


  AlphaS* mkAlphaS(const Info& info) {
    AlphaS* as = 0;
    const string itype = to_lower_copy(info.metadata("AlphaS_Type"));
    if (itype == "analytic") as = new AlphaS_Analytic();
    else if (itype == "ode") as = new AlphaS_ODE();
    else if (itype == "ipol") as = new AlphaS_Ipol();
    else throw FactoryError("Undeclared AlphaS requested: " + itype);
    // Configure the QCD params on this AlphaS
    if (info.has_key("AlphaS_OrderQCD")) as->setQCDorder(info.metadata_as<int>("AlphaS_OrderQCD"));
    if (!as->type() == "ipol") {
      if (info.has_key("MUp")) as->setQmass(1, info.metadata_as<double>("MUp"));
      if (info.has_key("MDown")) as->setQmass(2, info.metadata_as<double>("MDown"));
      if (info.has_key("MStrange")) as->setQmass(3, info.metadata_as<double>("MStrange"));
      if (info.has_key("MCharm")) as->setQmass(4, info.metadata_as<double>("MCharm"));
      if (info.has_key("MBottom")) as->setQmass(5, info.metadata_as<double>("MBottom"));
      if (info.has_key("MTop")) as->setQmass(6, info.metadata_as<double>("MTop"));
    }
    if (as->type() == "ode") {
      if( !(info.has_key("AlphaS_MZ")) || !(info.has_key("MZ")) )throw Exception("Requested ODE AlphaS without giving parameters for solving ODE.");
      as->setAlphaSMZ(info.metadata_as<double>("AlphaS_MZ"));
      as->setMZ(info.metadata_as<double>("MZ"));
    }
    if (as->type() == "analytic") {
      if( !(info.has_key("Lambda5")) && !(info.has_key("Lambda4")) &&
          !(info.has_key("Lambda3")) )throw Exception("Requested analytic AlphaS without setting a lambdaQCD.");
      if (info.has_key("Lambda3")) as->setLambda(3, info.metadata_as<double>("Lambda3"));
      if (info.has_key("Lambda4")) as->setLambda(4, info.metadata_as<double>("Lambda4"));
      if (info.has_key("Lambda5")) as->setLambda(5, info.metadata_as<double>("Lambda5"));
    }
    if (as->type() == "ipol") { /* populate interpolation vector */ }
    return as;
  }


}
