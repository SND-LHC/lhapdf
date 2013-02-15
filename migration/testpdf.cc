// Program to test LHAPDF6 PDF behaviour by writing out their values at lots of x and Q points

#include "LHAPDF/PDF.h"
#include "LHAPDF/Factories.h"
#include <iostream>
#include <fstream>
using namespace LHAPDF;
using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 2) {
    cerr << "You must specify a PDF set and member number" << endl;
    return 1;
  }

  const string setname = argv[1];
  const string smem = argv[2];
  const int imem = boost::lexical_cast<int>(smem);

  const PDF* pdf = mkPDF(setname, imem);
  vector<int> pids = pdf->flavors();

  foreach (int pid, pids) {
    const string spid = boost::lexical_cast<string>(pid);
    const string filename = setname + "_" + smem + "_" + spid + ".dat";
    ofstream f(filename.c_str());
    const double dx = 0.1; // 0.01
    const double dq2 = 0.1; // 0.01
    for (double log10x = -8; log10x < -0.1; log10x += dx) {
      for (double log10q2 = 1; log10q2 < 8; log10q2 += dq2) {
        const double x = pow(10, log10x);
        const double q2 = pow(10, log10q2);
        f << log10x << " " << log10q2 << " " << pdf->xfxQ2(pid, x, q2) << endl;
        // f << x << " " << q2 << " " << pdf->xfxQ2(pid, x, q2) << endl;
      }
    }
    f.close();
  }

  delete pdf;
  return 0;
}
