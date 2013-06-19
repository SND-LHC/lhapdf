// Example program to test alpha_s calculation

#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Factories.h"

#include <iostream>
#include <fstream>

using namespace LHAPDF;
using namespace std;


int main() {
  AlphaS_Analytic as_ana;

  as_ana.setQmass(1, 0.0017);
  as_ana.setQmass(2, 0.0041);
  as_ana.setQmass(3, 0.1);
  as_ana.setQmass(4, 1.29);
  as_ana.setQmass(5, 4.1);
  as_ana.setQmass(6, 172.5);

  as_ana.setLambda(3, 0.339);
  as_ana.setLambda(4, 0.296);
  as_ana.setLambda(5, 0.213);

  ofstream fa("alphas_ana.dat");//, fo("alphas_ode.dat");
  for (double log10q = 0; log10q < 2.5; log10q += 0.05) {
    const double q = pow(10, log10q);
    const double as_ana_q = as_ana.alphasQ(q);
    //cout << "alpha_s(Q=" << q << ")_ana = " << as_ana_q << endl;
    fa << q << " " << as_ana_q << endl;
    cout << "Q = " << setprecision(1) << fixed << q
         << ",  nf = " << as_ana.nf_Q(q) << endl;
    //const double as_ode_q = as_ode.alphasQ(q);
    //cout << "alpha_s(Q=" << q << ")_ODE = " << as_ode_q << endl;
    //fo << q << " " << as_ode_q << endl;
  }
  fa.close(); //fo.close();

  return 0;
}
