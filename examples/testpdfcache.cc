#include "LHAPDF/LHAPDF.h"
#include <memory>
using namespace std;

int main() {

  // Get a PDF set object and use it to get a vector of heap-allocated PDFs
  LHAPDF::PDFSet set("CT10nlo");
  const auto pdfs = set.mkPDFs<shared_ptr<const LHAPDF::PDF>>();
  for (double q = 10; q < 1e4; q *= 9) {
    for (double x = 1e-7; x < 1; x *= 10) {
      for (auto p : pdfs) {
        const double xf_g = p->xfxQ(21, 1e-3, 126.0);
        cout << q << ", " << x << ": " << xf_g << endl;
      }
    }
  }

  return 0;
}
