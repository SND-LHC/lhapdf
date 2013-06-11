// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"

namespace LHAPDF {


  /// Constructor from a set name and member ID.
  PDFInfo::PDFInfo(const std::string& setname, int member) {
    _setname = setname; _member = member;
    path searchpath = findFile(pdfmempath(setname, member));
    load(searchpath.string());
  }

  /// Constructor from an LHAPDF ID code.
  PDFInfo::PDFInfo(int lhaid) {
    const pair<string,int> setname_memid = lookupPDF(lhaid);
    if (setname_memid.second == -1)
      throw IndexError("Can't find a PDF with LHAPDF ID = " + to_str(lhaid));
    _setname = setname_memid.first; _member = setname_memid.second;
    path searchpath = pdfmempath(setname_memid.first, setname_memid.second);
    load(searchpath.string());
  }


  bool PDFInfo::has_key(const string& key) const {
    // cout << key << " in PDF: " << boolalpha << has_key_local(key) << endl;
    // cout << key << " in Set: " << boolalpha << getPDFSet(_setname).has_key(key) << endl;
    // cout << key << " in Cfg: " << boolalpha << getConfig().has_key(key) << endl;
    return has_key_local(key) || getPDFSet(_setname).has_key(key);
  }


  const std::string& PDFInfo::metadata(const string& key) const {
    if (has_key_local(key)) return metadata_local(key); //< value is defined locally
    return getPDFSet(_setname).metadata(key); //< fall back to the set-level info... or beyond
  }


  // /// @todo Only support loading via PDF set name and member ID, not explicit paths
  // /// @todo Replace the loading of the set metadata into the member info with set-level Info singletons
  // void Info::loadFull(const path& mempath) {
  //   // Extract the set name from the member data file path
  //   const path memberdata = findFile(mempath);
  //   if (memberdata.empty() || !exists(memberdata)) throw ReadError("Could not find PDF data file '" + mempath.string() + "'");
  //   const string memname = memberdata.filename().string(); //< Can use this to alternatively work out the set name...
  //   const path setdir = memberdata.parent_path();
  //   const string setname = setdir.filename().string();
  //   path setinfo = findpdfsetinfopath(setname);
  //   // Load the set info
  //   if (exists(setinfo)) load(setinfo.string());
  //   // Load the member info (possibly overriding the set-level metadata)
  //   load(memberdata.string());
  // }


}
