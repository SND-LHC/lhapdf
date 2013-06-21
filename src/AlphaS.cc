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
    /// Only mz and alphas_mz set to reasonable values...
    _lambda3 = -1;
    _lambda4 = -1;
    _lambda5 = -1;
    qcdorder = 0;
    mz = 91.19;
    alphas_mz = 0.118;
    _masstype = 0;
  }


  /// Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::nf_Q2(double q2) const {
    int nf = _nfmin;
    for (int it = _nfmin; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return (nf > _nfmax) ? _nfmax : nf;
  }


  /// Calculate a beta function given the number of active flavours
  double AlphaS::_beta(int i, int nf) const {
    if (i == 0) return (33 - 2*nf)/(12*M_PI);
    if (i == 1) return (153 - 19*nf)/(24*M_PI*M_PI);
    if (i == 2) return (2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf))/(128*M_PI*M_PI*M_PI);
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  /// Calculate beta functions given the number of active flavours
  vector<double> AlphaS::_betas(int nf) const {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 3; ++i) rtn.push_back(_beta(i, nf));
    return rtn;
  }

  /// Return the correct lambda for a given number of active flavours
  double AlphaS::_lambdaQCD(int nf) const {
    double lambda = -1;
    if(nf <= 3 && nf >= 0)lambda = _lambda3;
    else if(nf == 4)lambda = _lambda4;
    else if(nf > 4)lambda = _lambda5;
    if(lambda < 0)throw Exception("Invalid nf " + to_str(nf) + " requested for lambdaQCD");
    else return lambda;
  }

  /// Recalculate nfmax and nfmin after a new lambda has been set
  void AlphaS::_setFlavors() {
    if(_lambda3 < 0 && _lambda4 < 0 && _lambda5 < 0)throw Exception("No lambdas set");
    if(_lambda3 > 0 && _lambda4 < 0 && _lambda5 > 0)throw Exception("lambda3/5 set but not lambda4");
    _nfmin = (_lambda3 < 0) ? ( (_lambda4 < 0) ? 5 : 4) : 3;
    _nfmax = (_lambda5 < 0) ? ( (_lambda4 < 0) ? 3 : 4) : 5;
  }

  /// Set lambda_i && recalculate nfmax and nfmin
  void AlphaS::setLambda(int i, double lambda) {
    if(i != 3 && i != 4 && i != 5)throw Exception("Only lambda3/4/5 supported");
    if(i == 3)_lambda3 = lambda;
    if(i == 4)_lambda4 = lambda;
    if(i == 5)_lambda5 = lambda;
    _setFlavors();
  }

  /// Set a quark mass, and sort the mass vector accordingly
  void AlphaS::setQmass(double value) {
    if(_masstype == 2)throw Exception("Not allowed to mix the two ways of adding quark masses");
    _qmasses.push_back(value);
    sort( _qmasses.begin(), _qmasses.end() );
    _masstype = 1;
  }

  /// Set a quark mass, explicitly giving its id
  void AlphaS::setQmass(int id, double value) {
    if(_masstype == 1)throw Exception ("Not allowed to mix the two ways of adding quark masses");
    if(abs(id) > 6)throw Exception("Invalid id for quark given (should be 1-6).");
    if(_qmasses.size() != 6)_qmasses.resize(6);
    _qmasses[abs(id)-1] = value;
    _masstype = 2;
  }
}
