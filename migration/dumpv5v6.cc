/// PDF value dumping program, for both LHAPDF version 5 and 6
///
/// The version is detected automatically from the library against which the
/// executable is being compiled, so to build just do this:
/// g++ dumpv5v6.cc `lhapdf-config --cppflags -cxxflags --ldflags --libs` -o dumpv5orv6

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include <LHAPDF/LHAPDF.h>

#ifdef LHAPDF_MAJOR_VERSION
#define LHAVERSION LHAPDF_MAJOR_VERSION
#define XF(X, Q) pdf.xfxQ(flavor, X, Q)
#define AS(Q) pdf.alphasQ(Q)
#else
#define LHAVERSION 5
#define XF(X, Q) LHAPDF::xfx(X, Q, flavor)
#define AS(Q) LHAPDF::alphasPDF(Q)
#endif

#define _STR(s) #s
#define _VERSION(V) "v" _STR(V)
#define VERSION _VERSION(LHAVERSION)

const double MINLOGX = -10;
const double MINLOGQ = log10(1);
const double MAXLOGQ = log10(10000);
const double DX = 0.01;
const double DQ = 0.01;

int main(int argc, const char* argv[]) {

  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <setname> <setmember>" << std::endl;
    return 1;
  }
  const char* setname = argv[1];
  const int member = atoi(argv[2]);

  #if LHAVERSION > 5
  const LHAPDF::PDF& pdf = *LHAPDF::mkPDF(setname, member);
  #else
  string setfile = setname;
  if (setfile.find(".LH" == string::npos) setfile += ".LHgrid";
  LHAPDF::initPDFSetByName(setfile);
  LHAPDF::initPDF(member);
  #endif

  // Dump out points in (x,Q)
  for (int flavor = -6; flavor <= 6; ++flavor) {

    // x sampling for fixed Q
    double qs[] = {10, 50, 100, 200, 500, 1000, 2000, 5000};
    for (size_t iq = 0; iq < 8; ++iq) {
      const double q = qs[iq];
      std::stringstream filename;
      filename << "xf_" << VERSION << "_scanx_q" << q << "_" << flavor << ".dat";
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
      filename << "xf_" << VERSION << "_scanq_x" << x << "_" << flavor << ".dat";
      std::ofstream output(filename.str().c_str());
      for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
        const double q = pow(10, logQ);
        output << x << " " << q << " " << XF(x, q) << std::endl;
      }
      output.close();
    }

    // alpha_s sampling in Q
    std::stringstream as_filename;
    as_filename << "as_" << VERSION << ".dat";
    std::ofstream as_output(as_filename.str().c_str());
    for (double logQ = MINLOGQ; logQ <= MAXLOGQ; logQ += DQ) {
      const double q = pow(10, logQ);
      as_output << q << " " << AS(q) << std::endl;
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
