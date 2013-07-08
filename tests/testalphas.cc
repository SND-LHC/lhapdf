 // Example program to test alpha_s calculation

#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>

using namespace LHAPDF;
using namespace std;


int main() {
  AlphaS_Analytic as_ana;
  as_ana.setOrderQCD(4);
  as_ana.setQMass(1, 0.0017);
  as_ana.setQMass(2, 0.0041);
  as_ana.setQMass(3, 0.1);
  as_ana.setQMass(4, 1.29);
  as_ana.setQMass(5, 4.1);
  as_ana.setQMass(6, 172.5);
  // as_ana.setLambda(2, 0.400);
  as_ana.setLambda(3, 0.339);
  as_ana.setLambda(4, 0.296);
  as_ana.setLambda(5, 0.213);

  AlphaS_ODE as_ode;
  as_ode.setMZ(91);
  as_ode.setAlphaSMZ(0.118);
  as_ode.setQMass(1, 0.0017);
  as_ode.setQMass(2, 0.0041);
  as_ode.setQMass(3, 0.1);
  as_ode.setQMass(4, 1.29);
  as_ode.setQMass(5, 4.1);
  as_ode.setQMass(6, 172.5);

  AlphaS_Ipol as_ipol;
  vector<double> qs;
  qs.push_back(1);
  qs.push_back(91);
  vector<double> alphas;
  alphas.push_back(0.5);
  alphas.push_back(0.118);
  as_ipol.setQValues(qs);
  as_ipol.setAlphaSValues(alphas);

  PDF* pdf = mkPDF("CT10nlo", 0);

  ofstream fa("alphas_ana.dat"), fo("alphas_ode.dat"), fi("alphas_ipol.dat"), fc("alphas_ct10nlo.dat");
  for (double log10q = -0.5; log10q < 3; log10q += 0.05) {
    const double q = pow(10, log10q);
    const double as_ana_q = as_ana.alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_ana = " << as_ana_q << endl;
    fa << q << " " << as_ana_q << endl;
    const double as_ode_q = as_ode.alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_ODE = " << as_ode_q << endl;
    fo << q << " " << as_ode_q << endl;
    const double as_ipol_q = as_ipol.alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_ipol = " << as_ipol_q << endl;
    fi << q << " " << as_ipol_q << endl;
    const double as_ct10_q = pdf->alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_CT10nlo = " << as_ct10_q << endl;
    fc << q << " " << as_ct10_q << endl;
  }
  fa.close(); fo.close(); fi.close(); fc.close();

  delete pdf;
  return 0;
}
