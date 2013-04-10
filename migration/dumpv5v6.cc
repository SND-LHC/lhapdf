/// PDF value dumping program, for both LHAPDF version 5 and 6
///
/// The version is specified via the LHAVERSION macro definition, e.g.
///   g++ dumpv5v6.cc -DLHAVERSION=5 -o dumpv5 `lhapdf-config --cppflags -cxxflags --ldflags --libs`
///   g++ dumpv5v6.cc -DLHAVERSION=6 -o dumpv6 `lhapdf-config --cppflags -cxxflags --ldflags --libs`

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include <LHAPDF/LHAPDF.h>

#ifdef LHAPDF_MAJOR_VERSION
#define LHAVERSION LHAPDF_MAJOR_VERSION
#define XF(X, Q) pdf.xfxQ(flavor, X, Q)
#else
#define LHAVERSION 5
#define XF(X, Q) LHAPDF::xfx(X, Q, id)
#endif

#define _STR(s) #s
#define _OUTPUT_PREFIX(V) "xf_v" _STR(V)
#define OUTPUT_PREFIX _OUTPUT_PREFIX(LHAVERSION)

const double MINLOGX = -10;
const double MINLOGQ = log10(10);
const double MAXLOGQ = log10(5000);
const double DX = 0.01;
const double DQ = 0.01;

int main(int argc, const char* argv[]) {

  if (argc < 3) {
    std::cerr << argv[0] << " <setname> <setmember>" << std::endl;
    return 1;
  }
  const char* setname = argv[1];
  const int member = atoi(argv[2]);

  #if LHAVERSION > 5
  const LHAPDF::PDF& pdf = *LHAPDF::mkPDF(setname, member);
  #else
  LHAPDF::initPDFSetByName(setname);
  LHAPDF::initPDF(member);
  #endif

  // Dump out points in (x,Q)
  int flavors[] = {-5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5};
  for (int i = 0; i < 11; ++i) {
    const int flavor = flavors[i];
    #if LHAVERSION > 5
    const int id = i;
    #else
    const int id = i-5;
    #endif

    // x sampling for fixed Q
    double qs[] = {10, 50, 100, 200, 500, 1000, 2000, 5000};
    for (size_t iq = 0; iq < 8; ++iq) {
      const double q = qs[iq];
      std::stringstream filename;
      filename << OUTPUT_PREFIX << "_scanx_q" << q << "_" << flavor << ".dat";
      std::ofstream output(filename.str().c_str());
      for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
        const double x  = pow(10, logX);
        output << x << " " << q << " " << XF(x, q) << std::endl;
      }
      output.close();
    }

    // Q sampling for fixed x
    double xs[] = {1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.2, 0.5, 0.8};
    for (size_t ix = 0; ix < 8; ++ix) {
      const double x = xs[ix];
      std::stringstream filename;
      filename << OUTPUT_PREFIX << "_scanq_x" << x << "_" << flavor << ".dat";
      std::ofstream output(filename.str().c_str());
      for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
        const double q = pow(10, logQ);
        output << x << " " << q << " " << XF(x, q) << std::endl;
      }
      output.close();
    }

    // // (x,Q) 2D sampling
    // std::stringstream filename;
    // filename << OUTPUT_PREFIX << "_scanxq" << flavor << ".dat";
    // std::ofstream output(filename.str().c_str());
    // for (double logX = MINLOGX; logX <= 0.0; logX += DX) {
    //   for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
    //     const double x  = pow(10, logX);
    //     const double q = pow(10, logQ);
    //     output << x << " " << q << " " << XF(x, q) << std::endl;
    //   }
    // }
    // output.close();

  }

  return 0;
}
