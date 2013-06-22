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
    // Only mz and alphas_mz set to reasonable values...
    _lambdas.clear();
    qcdorder = 0;
    mz = 91.19;
    alphas_mz = 0.118;
    _masstype = 0;
  }


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::nf_Q2(double q2) const {
    int nf = _nfmin;
    for (int it = _nfmin; it < (int)_qmasses.size(); ++it)
      if (q2 > sqr(_qmass(it))) nf = it;
    return (nf > _nfmax) ? _nfmax : nf;
  }


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

  // Return the correct lambda for a given number of active flavours
  double AlphaS::_lambdaQCD(int nf) const {
    double lambda = _lambdas[nf-1];
    if(lambda == 0)throw Exception("Invalid nf " + to_str(nf) + " requested for lambdaQCD");
    return lambda;
  }

  // Recalculate nfmax and nfmin after a new lambda has been set
  // Relies on vector<double> initialising all elements to 0 by
  // default, I *think* this should work even though I am comparing
  // floats since 0 has a precise representation
  void AlphaS::_setFlavors() {
    for(unsigned int it = 0; it < _lambdas.size(); ++it) {
      if(_lambdas[it] != 0.){_nfmin = it+1; break;}
    }
    for(unsigned int it = _lambdas.size() - 1; it >= 0; --it) {
      if(_lambdas[it] != 0.){_nfmax = it+1; break;}
    }
    if(_nfmin != _nfmax) {
      for(int it = _nfmin; it < _nfmax - 1; ++it){
        if(_lambdas[it] == 0)throw Exception ("Need to set intermediate lambdas.");
      }
    }
  }

  // Set lambda_i && recalculate nfmax and nfmin
  void AlphaS::setLambda(unsigned int i, double lambda) {
    if(_lambdas.size() < i)_lambdas.resize(i);
    _lambdas[i-1] = lambda;
    _setFlavors();
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
