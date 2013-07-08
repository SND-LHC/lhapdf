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


  Info* mkPDFInfo(const std::string& setname, int member) {
    return new Info(findpdfmempath(setname, member));
  }

  Info* mkPDFInfo(int lhaid) {
    const pair<string,int> setname_memid = lookupPDF(lhaid);
    return mkPDFInfo(setname_memid.first, setname_memid.second);
  }


  PDF* mkPDF(const string& setname, int member) {
    // First create an Info object to work out what format of PDF this is:
    Info info(findpdfmempath(setname, member));
    const string fmt = info.get_entry("Format");
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

    const string itype = to_lower_copy(info.get_entry("AlphaS_Type"));
    if (itype == "analytic") as = new AlphaS_Analytic();
    else if (itype == "ode") as = new AlphaS_ODE();
    else if (itype == "ipol") as = new AlphaS_Ipol();
    else throw FactoryError("Undeclared AlphaS requested: " + itype);

    // Configure the QCD params on this AlphaS
    if (info.has_key("AlphaS_OrderQCD")) as->setOrderQCD(info.get_entry_as<int>("AlphaS_OrderQCD"));
    if (info.has_key("MUp")) as->setQMass(1, info.get_entry_as<double>("MUp"));
    if (info.has_key("MDown")) as->setQMass(2, info.get_entry_as<double>("MDown"));
    if (info.has_key("MStrange")) as->setQMass(3, info.get_entry_as<double>("MStrange"));
    if (info.has_key("MCharm")) as->setQMass(4, info.get_entry_as<double>("MCharm"));
    if (info.has_key("MBottom")) as->setQMass(5, info.get_entry_as<double>("MBottom"));
    if (info.has_key("MTop")) as->setQMass(6, info.get_entry_as<double>("MTop"));

    // Required parameter settings for each calculation mode
    if (as->type() == "ode") {
      /// @todo Handle FFNS / VFNS
      if (!info.has_key("AlphaS_MZ") || !info.has_key("MZ") )
        throw MetadataError("Requested ODE AlphaS but the required parameters are not defined.");
      as->setAlphaSMZ(info.get_entry_as<double>("AlphaS_MZ"));
      as->setMZ(info.get_entry_as<double>("MZ"));
      if (info.has_key("AlphaS_Qs")) {
        AlphaS_ODE* as_o = dynamic_cast<AlphaS_ODE*>(as);
        if (info.has_key("AlphaS_Qs")) as_o->setQValues( info.get_entry_as< vector<double> >("AlphaS_Qs"));
      }
    }
    else if (as->type() == "analytic") {
      /// @todo Handle FFNS / VFNS
      if (!info.has_key("Lambda5") && !info.has_key("Lambda4") && !info.has_key("Lambda3") )
        throw MetadataError("Requested analytic AlphaS but the required parameters are not defined.");
      if (info.has_key("Lambda3")) as->setLambda(3, info.get_entry_as<double>("Lambda3"));
      if (info.has_key("Lambda4")) as->setLambda(4, info.get_entry_as<double>("Lambda4"));
      if (info.has_key("Lambda5")) as->setLambda(5, info.get_entry_as<double>("Lambda5"));
    }
    else if (as->type() == "ipol") {
      if (!info.has_key("AlphaS_Qs") || !info.has_key("AlphaS_Vals") )
        throw MetadataError("Requested ipol AlphaS but the required parameters are not defined.");
      AlphaS_Ipol* as_i = dynamic_cast<AlphaS_Ipol*>(as);
      if (info.has_key("AlphaS_Qs")) as_i->setQValues( info.get_entry_as< vector<double> >("AlphaS_Qs"));
      if (info.has_key("AlphaS_Vals")) as_i->setAlphaSValues( info.get_entry_as< vector<double> >("AlphaS_Vals"));
    }

    return as;
  }


}
