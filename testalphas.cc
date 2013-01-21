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
  for (double log10q = 0.2; log10q < 3; log10q += 0.05) {
    const double q = pow(10, log10q);
    cout << "Q = " << q << ", alpha_s(Q)_ana = " << as_ana.alphasQ(q) << endl;
    fa << q << " " << as_ana.alphasQ(q) << endl;
    // cout << "Q = " << q << ", alpha_s(Q)_ODE = " << as_ode.alphasQ(q) << endl;
    // fo << q << " " << as_ode.alphasQ(q) << endl;
  }
  fa.close(); fo.close();

  return 0;
}
