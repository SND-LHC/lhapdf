#include "LHAPDF/LHAPDF.h"

#include <string>
#include <vector>
#include <map>
using namespace std;

int main() {
  LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10nlo", 0);
  double xf_g = pdf->xfxQ(21, 1e-3, 126.0);
  cout << xf_g << endl;
  map<int, double> xfs = pdf->xfxQ(1e-3, 126.0);

  size_t num_mems = pdf->numMembers();
  delete pdf;

  typedef unique_ptr<LHAPDF::PDF> PdfPtr;
  vector<PdfPtr> pdfs;
  for (size_t i = 0; i < num_mems; ++i)
    pdfs.push_back( PdfPtr(LHAPDF::mkPDF("CT10nlo", i)) );
  for (const auto& p : pdfs) {
    double xf_g = p->xfxQ(21, 1e-3, 126.0);
    cout << xf_g << endl;
  }

  return 0;
}
