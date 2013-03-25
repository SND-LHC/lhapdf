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
    int flavor = flavors[id];

    std::stringstream filename;
    filename << "xf_v6_" << flavor << ".dat";
    std::ofstream output( filename.str().c_str() );

    // x sampling for fixed Q
    const double q = 500;
    for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
      double x  = pow(10, logX);
      output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
    }

    // // Q sampling for fixed x
    // const double x = 1e-3;
    // for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
    //   double q = pow(10, logQ);
    //   output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
    // }

    // // (x,Q) 2D sampling
    // for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
    //   for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
    //     const double x  = pow(10, logX);
    //     const double q = pow(10, logQ);
    //     output << x << " " << q << " " << pdf.xfxQ(flavor, x, q) << std::endl;
    //   }
    // }

    output.close();
  }

  return 0;
}
