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
    /// @todo Should we instead use explicitly invalid defaults?
    setQmass(1, -1);
    setQmass(2, -1);
    setQmass(3, -1);
    setQmass(4, -1);
    setQmass(5, -1);
    setQmass(6, -1);
    _lambda3 = -1;
    _lambda4 = -1;
    _lambda5 = -1;
    qcdorder = 0;
    mz = 91.19;
    alphas_mz = 0.118;
  }


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::nf_Q2(double q2) const {
    int nf = _nfmin;
    for (int it = _nfmin; it <= 5; ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return (nf > _nfmax) ? _nfmax : nf;
  }


  // Calculate a beta function given the number of active flavours
  double AlphaS::_beta(int i, int nf) const {
    if (i == 0) return (33 - 2*nf)/(12*M_PI);
    if (i == 1) return (153 - 19*nf)/(24*M_PI*M_PI);
    if (i == 2) return (2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf))/(128*M_PI*M_PI*M_PI);
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::_betas(int nf) const {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 3; ++i) rtn.push_back(_beta(i, nf));
    return rtn;
  }

  // Return the correct lambda for a given number of active flavours
  double AlphaS::_lambdaQCD(int nf) const {
    double lambda = -1;
    if(nf <= 3 && nf >= 0)lambda = _lambda3;
    else if(nf == 4)lambda = _lambda4;
    else if(nf > 4)lambda = _lambda5;
    if(lambda < 0)throw Exception("Invalid nf " + to_str(nf) + " requested for lambdaQCD");
    else return lambda;
  }

  void AlphaS::_setFlavors() {
    if(_lambda3 < 0 && _lambda4 < 0 && _lambda5 < 0)throw Exception("No lambdas set");
    if(_lambda3 > 0 && _lambda4 < 0 && _lambda5 > 0)throw Exception("lambda3/5 set but not lambda4");
    _nfmin = (_lambda3 < 0) ? ( (_lambda4 < 0) ? 5 : 4) : 3;
    _nfmax = (_lambda5 < 0) ? ( (_lambda4 < 0) ? 3 : 4) : 5;
  }

  void AlphaS::setLambda(int i, double lambda) {
    if(i != 3 && i != 4 && i != 5)throw Exception("Only lambda3/4/5 supported");
    if(i == 3)_lambda3 = lambda;
    if(i == 4)_lambda4 = lambda;
    if(i == 5)_lambda5 = lambda;
    _setFlavors();
  }

}
