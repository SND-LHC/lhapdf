// Example program to test PDF grid format reading and interpolation

#include "LHAPDF/AlphaS.h"
#include <iostream>
#include <fstream>
using namespace LHAPDF;
using namespace std;

int main() {
  const AlphaS_ODE as_ode;
  const AlphaS_Analytic as_ana;

  ofstream fa("alphas_ana.dat"), fo("alphas_ode.dat");
  for (double log10q = 0.8; log10q < 2.5; log10q += 0.05) {
    const double q = pow(10, log10q);
    const double as_ana_q = as_ana.alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_ana = " << as_ana_q << endl;
    fa << q << " " << as_ana_q << endl;
    const double as_ode_q = as_ode.alphasQ(q);
    cout << "alpha_s(Q=" << q << ")_ODE = " << as_ode_q << endl;
    fo << q << " " << as_ode_q << endl;
  }
  fa.close(); fo.close();

  return 0;
}
