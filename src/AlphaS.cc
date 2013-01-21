#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  // Base class constructor for default param setup
  AlphaS::AlphaS() {
    _qmasses[0] = 0.0017;
    _qmasses[1] = 0.0041;
    _qmasses[2] = 0.1;
    _qmasses[3] = 1.29;
    _qmasses[4] = 4.1;
    _qmasses[5] = 172.5;
    _mz = 91.1876;
    _alphas_mz = 0.118;
    /// @todo Use better, or more explicitly invalid, defaults
    _lambda4 = 0.215;
    _lambda5 = 0.165;
    _order = 1;
  }


  // Calculate the number of active quark flavours at energy scale Q2
  int AlphaS::nf_Q2(double q2) const {
    for (int nf = 2; nf <= 4; ++nf)
      if (q2 < sqr(_qmasses[nf])) return nf;
    return 5;
  }


  // Calculate a beta function given the number of active flavours
  double AlphaS::beta(int i, int nf) const {
    if (i == 0) return 11 - (2 / 3.0)*nf;
    if (i == 1) return 51 - (19 / 3.0)*nf;
    if (i == 2) return 2857 - (5033 / 9.0)*nf + (325 / 27.0)*sqr(nf);
    throw Exception("Invalid index " + to_str(i) + " for requested beta function");
  }


  // Calculate beta functions given the number of active flavours
  vector<double> AlphaS::betas(int nf) const {
    vector<double> rtn; rtn.reserve(3);
    for (int i = 0; i < 3; ++i) rtn.push_back(beta(i, nf));
    return rtn;
  }


}
