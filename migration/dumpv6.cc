#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include <LHAPDF/PDF.h>
#include <LHAPDF/Factories.h>

const double MINLOGX = -10;
const double MINLOGQ = log10(10);
const double MAXLOGQ = log10(5000);
const double DX = 0.1;
const double DQ = 0.1;

int main(int argc, const char* argv[]) {

  if (argc < 3) {
    std::cerr << argv[0] << " <setname> <setmember>" << std::endl;
    return 1;
  }

  const char* setname = argv[1];
  const int setmember = atoi(argv[2]);
  const LHAPDF::PDF& pdf = *LHAPDF::mkPDF(setname, setmember);

  // Dump out points in (x,Q)
  int flavors[] = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
  for (int id = 0; id < 11; ++id) {
    const int flavor = flavors[id];

    // x sampling for fixed Q
    double qs[] = {10, 50, 100, 200, 500, 1000, 2000, 5000};
    for (size_t iq = 0; iq < 8; ++iq) {
      const double q = qs[iq];
      std::stringstream filename;
      filename << "xf_v6_scanQ" << q << "_" << flavor << ".dat";
      std::ofstream output(filename.str().c_str());
      for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
        const double x  = pow(10, logX);
        output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
      }
      output.close();
    }

    // Q sampling for fixed x
    double xs[] = {1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.2, 0.5, 0.8};
    for (size_t ix = 0; ix < 8; ++ix) {
      const double x = xs[ix];
      std::stringstream filename;
      filename << "xf_v6_scanx" << x << "_" << flavor << ".dat";
      std::ofstream output(filename.str().c_str());
      for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
        const double q = pow(10, logQ);
        output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
      }
      output.close();
    }

    // (x,Q) 2D sampling
    std::stringstream filename;
    filename << "xf_v6_scanxQ" << flavor << ".dat";
    std::ofstream output(filename.str().c_str());
    for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
      for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
        const double x  = pow(10, logX);
        const double q = pow(10, logQ);
        output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
      }
    }
    output.close();
  }

  return 0;
}
