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
    setQmass(1, 0.0017);
    setQmass(2, 0.0041);
    setQmass(3, 0.1);
    setQmass(4, 1.29);
    setQmass(5, 4.1);
    setQmass(6, 172.5);
    mz = 91.1876;
    alphas_mz = 0.118;
    lambda5 = 0.213;
    lambda4 = 0.296;
    lambda3 = 0.339;
    qcdorder = 1;
    maxnf = 6;
  }


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::nf_Q2(double q2) const {
    for (int nf = 2; nf <= 4; ++nf)
      if (q2 < sqr(qmass(nf+1))) return nf;
    return 5;
  }


  // Calculate a beta function given the number of active flavours
  double AlphaS::beta(int i, int nf) const {
    if (i == 0) return (33 - 2*nf)/(12*M_PI);
    if (i == 1) return (153 - 19*nf)/(24*M_PI*M_PI);
    if (i == 2) return (2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf))/(128*M_PI*M_PI*M_PI);
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::betas(int nf) const {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 3; ++i) rtn.push_back(beta(i, nf));
    return rtn;
  }


}
