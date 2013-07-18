// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Base class constructor for default param setup
  AlphaS::AlphaS() {
    _qcdorder = 4;
  }


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::numFlavorsQ2(double q2) const {
    /// @todo Respect the NumFlavors + FlavorScheme metadata
    int nf = 0;
    for (int it = 0; it < (int)_qmasses.size(); ++it)
      // if (q2 > sqr(_qmasses[it]) && _qmasses[it] != 0) nf = it + 1;
      if (q2 > sqr(_qmasses[it])) nf = it + 1; else break;
    return nf;
  }


  // Calculate a beta function given the number of active flavours
  double AlphaS::_beta(int i, int nf) const {
    if (i == 0) return (double) 0.875352187 - 0.053051647*nf; //(33 - 2*nf)/(12*M_PI)
    if (i == 1) return (double) 0.6459225457 - 0.0802126037*nf; //(153 - 19*nf)/(24*sqr(M_PI))
    if (i == 2) return (double) 0.719864327 - 0.140904490*nf + 0.00303291339*nf*nf; //(2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf))/(128*sqr(M_PI)*M_PI)
    if (i == 3) return (double) 1.172686 - 0.2785458*nf + 0.01624467*nf*nf + 0.0000601247*nf*nf*nf;
    // ( (149753/6.) + 3564*ZETA_3 - ((1078361/162.) + (6502/27.)*ZETA_3)*nf +
    // ((50065/162.) + (6472/81.)*ZETA_3)*sqr(nf) + (1093/729.)*sqr(nf)*nf)/(256*sqr(M_PI)*sqr(M_PI))

    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::_betas(int nf) const {
    vector<double> rtn; rtn.reserve(4);
    for (int i = 0; i < 4; ++i) rtn.push_back(_beta(i, nf));
    return rtn;
  }


  // Set a quark mass, explicitly giving its ID
  void AlphaS::setQuarkMass(int id, double value) {
    if (abs(id) > 6 || id == 0)
      throw Exception("Invalid id for quark given (should be 1-6).");
    if (_qmasses.size() != 6) _qmasses.resize(6);
    _qmasses[abs(id)-1] = value;
  }


  // Get a quark mass by ID
  double AlphaS::quarkMass(int id) const {
    if (abs(id) > 6 || id == 0)
       throw Exception("Invalid id for quark given (should be 1-6).");
    return _qmasses[abs(id)-1];
  }


}
