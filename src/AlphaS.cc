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
    _qcdorder = 0;
    _masstype = 0;
  }

  AlphaS::~AlphaS() {}

  // Calculate a beta function given the number of active flavours
  double AlphaS::_beta(int i, int nf) const {
    if (i == 0) return (double) (33 - 2*nf)/(12*M_PI);
    if (i == 1) return (double) (153 - 19*nf)/(24*sqr(M_PI));
    if (i == 2) return (double) (2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf))/(128*sqr(M_PI)*M_PI);
    if (i == 3) return (double) ( (149753/6.) + 3564*ZETA_3 - ((1078361/162.) + (6502/27.)*ZETA_3)*nf +
                                ((50065/162.) + (6472/81.)*ZETA_3)*sqr(nf) + (1093/729.)*sqr(nf)*nf)/(256*sqr(nf)*sqr(nf));
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::_betas(int nf) const {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 4; ++i) rtn.push_back(_beta(i, nf));
    return rtn;
  }

  // Set a quark mass, and sort the mass vector accordingly
  void AlphaS::setQmass(double value) {
    if(_masstype == 2)throw Exception("Not allowed to mix the two ways of adding quark masses");
    _qmasses.push_back(value);
    sort( _qmasses.begin(), _qmasses.end() );
    _masstype = 1;
  }

  // Set a quark mass, explicitly giving its id
  void AlphaS::setQmass(int id, double value) {
    if(_masstype == 1)throw Exception ("Not allowed to mix the two ways of adding quark masses");
    if(abs(id) > 6 || id == 0)throw Exception("Invalid id for quark given (should be 1-6).");
    if(_qmasses.size() != 6)_qmasses.resize(6);
    _qmasses[abs(id)-1] = value;
    _masstype = 2;
  }
}
