#pragma once

#include "LHAPDF/Paths.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Exceptions.h"

namespace LHAPDF {


  /// Get the singleton LHAPDF set ID -> PDF index map
  std::map<int, std::string>& getPDFIndex() {
    static map<int, string> _lhaindex;
    if (_lhaindex.empty()) { // The map needs to be populated first
      path indexpath = findFile("pdfsets.index");
      if (indexpath.empty()) throw ReadError("Could not find a pdfsets.index file");
      try {
        ifstream file(indexpath.c_str());
        string line;
        while (getline(file, line)) {
          istringstream tokens(line);
          int id; string setname;
          tokens >> id;
          tokens >> setname;
          // cout << id << " -> " << _lhaindex[id] << endl;
          _lhaindex[id] = setname;
        }
      } catch (const std::exception& ex) {
        throw ReadError("Trouble when reading " + indexpath.native() + ": " + ex.what());
      }
    }
    return _lhaindex;
  }


  /// Look up a PDF set name and member ID by the LHAPDF ID code
  ///
  /// The set name and member ID are returned as an std::pair.
  /// If lookup fails, a pair ("", -1) is returned.
  inline pair<std::string, int> lookupPDF(int lhaid) {
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


  /// Look up the LHAPDF index from the set name and member ID.
  ///
  /// If lookup fails, -1 is returned, otherwise the LHAPDF ID code.
  /// NB. This function is relatively slow, since it requires std::map reverse lookup.
  inline int lookupSetID(const std::string& setname, int memid) {
    // const map<int, string>& = getPDFIndex();
    typedef pair<int, string> MapPair;
    foreach (const MapPair& id_name, getPDFIndex()) {
      if (id_name.second == setname) return id_name.first + memid;
    }
    return -1; //< failure value
  }


}
