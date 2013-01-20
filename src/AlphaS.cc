#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::_nf_Q2(double q2) {
    int nf = 2;
    for ( ; nf < 6; ++nf) {
      if (q2 < sqr(_qmasses[nf])) break; // check against *next* mass
    }
    return nf; // max value = 5, even for q > mt
  }


  // Calculate a beta function given the number of active flavours
  double AlphaS::beta(int i, int nf) {
    if (i == 0) return 11 - (2 / 3.0)*nf;
    if (i == 1) return 51 - (19 / 3.0)*nf;
    if (i == 2) return 2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf);
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::betas(int nf) {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 3; ++i) rtn.push_back(beta(i, nf));
    return rtn;
  }


}
