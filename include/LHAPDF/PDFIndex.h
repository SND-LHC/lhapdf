// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2016 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_PDFIndex_H
#define LHAPDF_PDFIndex_H

#include "LHAPDF/Paths.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Exceptions.h"

namespace LHAPDF {


  /// @name Functions for PDF lookup by LHAPDF ID index file
  //@{

  /// Get the singleton LHAPDF set ID -> PDF index map
  inline std::map<int, std::string>& getPDFIndex() {
    static map<int, string> _lhaindex;
    if (_lhaindex.empty()) { // The map needs to be populated first
      string indexpath = findFile("pdfsets.index");
      if (indexpath.empty()) throw ReadError("Could not find a pdfsets.index file");
      try {
        ifstream file(indexpath.c_str());
        string line;
        while (getline(file, line)) {
          line = trim(line);
          if (line.empty() || line.find("#") == 0) continue;
          istringstream tokens(line);
          int id; string setname;
          tokens >> id;
          tokens >> setname;
          // cout << id << " -> " << _lhaindex[id] << endl;
          _lhaindex[id] = setname;
        }
      } catch (const std::exception& ex) {
        throw ReadError("Trouble when reading " + indexpath + ": " + ex.what());
      }
    }
    return _lhaindex;
  }


  /// Look up a PDF set name and member ID by the LHAPDF ID code
  ///
  /// The set name and member ID are returned as an std::pair.
  /// If lookup fails, a pair ("", -1) is returned.
  inline std::pair<std::string, int> lookupPDF(int lhaid) {
    map<int, string>::iterator it = getPDFIndex().upper_bound(lhaid);
    string rtnname = "";
    int rtnmem = -1;
    if (it != getPDFIndex().begin()) {
      --it; // upper_bound (and lower_bound) return the entry *above* lhaid: we need to step back
      rtnname = it->second; // name of the set that contains this ID
      rtnmem = lhaid - it->first; // the member ID is the offset from the lookup ID
    }
    return make_pair(rtnname, rtnmem);
  }


  /// @brief Decode a single PDF member ID string into a setname,memid pair
  ///
  /// @note A trivial <SET,MEM> decoding rather than a "rea; lookup", for convenience & uniformity.
  inline std::pair<std::string,int> lookupPDF(const std::string& pdfstr) {
    int nmem = 0;
    const size_t slashpos = pdfstr.find("/");
    const string setname = trim(pdfstr.substr(0, slashpos));
    try {
      if (slashpos != string::npos) {
        const string smem = pdfstr.substr(slashpos+1);
        nmem = lexical_cast<int>(smem);
      }
    } catch (...) {
      throw UserError("Could not parse PDF identity string " + pdfstr);
    }
    return make_pair(setname, nmem);
  }


  /// Look up the member's LHAPDF index from the set name and member ID.
  ///
  /// If lookup fails, -1 is returned, otherwise the LHAPDF ID code.
  /// NB. This function is relatively slow, since it requires std::map reverse lookup.
  inline int lookupLHAPDFID(const std::string& setname, int nmem) {
    // const map<int, string>& = getPDFIndex();
    typedef pair<int, string> MapPair;
    for (const MapPair& id_name : getPDFIndex()) {
      if (id_name.second == setname) return id_name.first + nmem;
    }
    return -1; //< failure value
  }


  /// Look up the member's LHAPDF index from a setname/member string.
  inline int lookupLHAPDFID(const std::string& setname_nmem) {
    const pair<string,int> idpair = lookupPDF(setname_nmem);
    return lookupLHAPDFID(idpair.first, idpair.second);
  }

  //@}

}
#endif
