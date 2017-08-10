// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2016 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
using namespace std;

namespace LHAPDF {


  double PDF::xfxQ2(int id, double x, double q2) const {
    // Physical x range check
    if (!inPhysicalRangeX(x)) {
      throw RangeError("Unphysical x given: " + to_str(x));
    }
    // Physical Q2 range check
    if (!inPhysicalRangeQ2(q2)) {
      throw RangeError("Unphysical Q2 given: " + to_str(q2));
    }
    // Treat PID = 0 as always equivalent to a gluon: query as PID = 21
    const int id2 = (id != 0) ? id : 21; //< @note Treat 0 as an alias for 21
    // Undefined PIDs
    if (!hasFlavor(id2)) return 0.0;
    // Call the delegated method in the concrete PDF object to calculate the in-range value
    double xfx = _xfxQ2(id2, x, q2);
    // Apply positivity forcing at the enabled level
    switch (forcePositive()) {
    case 0: break;
    case 1: if (xfx < 0) xfx = 0; break;
    case 2: if (xfx < 1e-10) xfx = 1e-10; break;
    default: throw LogicError("ForcePositive value not in expected range!");
    }
    // Return
    return xfx;
  }


  void PDF::xfxQ2(double x, double q2, std::map<int, double>& rtn) const {
    rtn.clear();
    for (int id : flavors()) rtn[id] = xfxQ2(id, x, q2);
  }


  void PDF::xfxQ2(double x, double q2, std::vector<double>& rtn) const {
    rtn.clear();
    rtn.resize(13);
    for (int i = 0; i < 13; ++i) {
      const int id = i-6; // PID = 0 is automatically treated as PID = 21
      rtn[i] = xfxQ2(id, x, q2);
    }
  }


  std::map<int, double> PDF::xfxQ2(double x, double q2) const {
    std::map<int, double> rtn;
    xfxQ2(x, q2, rtn);
    return rtn;
  }


  void PDF::print(std::ostream& os, int verbosity) const {
    stringstream ss;
    if (verbosity > 0) {
      ss << set().name() << " PDF set, member #" << memberID()
         << ", version " << dataversion();
      if (lhapdfID() > 0)
        ss << "; LHAPDF ID = " << lhapdfID();
    }
    if (verbosity > 2 && set().description().size() > 0)
      ss << "\n" << set().description();
    if (verbosity > 1 && description().size() > 0)
      ss << "\n" << description();
    if (verbosity > 2)
      ss << "\n" << "Flavor content = " << to_str(flavors());
    os << ss.str() << endl;
  }


  int PDF::lhapdfID() const {
    //return set().lhapdfID() + memberID()
    /// @todo Add failure tolerance if pdfsets.index not found
    try {
      return lookupLHAPDFID(_setname(), memberID());
    } catch (const Exception&) {
      return -1; //< failure
    }
  }


  double PDF::quarkMass(int id) const {
    const unsigned int aid = std::abs(id);
    if (aid == 0 || aid > 6) return -1;
    const static string QNAMES[] = {"Down", "Up", "Strange", "Charm", "Bottom", "Top"}; ///< @todo Centralise?
    const size_t qid = aid - 1;
    const string qname = QNAMES[qid];
    return info().get_entry_as<double>("M" + qname, -1);
  }


  double PDF::quarkThreshold(int id) const {
    const unsigned int aid = std::abs(id);
    if (aid == 0 || aid > 6) return -1;
    const static string QNAMES[] = {"Down", "Up", "Strange", "Charm", "Bottom", "Top"}; ///< @todo Centralise?
    const size_t qid = aid - 1;
    const string qname = QNAMES[qid];
    return info().get_entry_as<double>("Threshold" + qname, quarkMass(id));
  }


}
