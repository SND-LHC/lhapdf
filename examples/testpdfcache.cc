#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/LogBicubicInterpolator.h"
#include <memory>
using namespace std;

int main() {

  // Get a PDF set object and use it to get a vector of heap-allocated PDFs
  LHAPDF::PDFSet set("CT10nlo");
  const auto pdfs = set.mkPDFs<shared_ptr<const LHAPDF::PDF>>();

  cout << "Original cache size = " << LHAPDF::LogBicubicInterpolator::XCaches::N << endl;
  cout << pdfs[0]->xfxQ(21, 1e-3, 126.0) << endl;

  LHAPDF::LogBicubicInterpolator::XCaches::setup(16);
  cout << "New cache size = " << LHAPDF::LogBicubicInterpolator::XCaches::N << endl;
  cout << pdfs[0]->xfxQ(21, 1e-3, 126.0) << endl;

  for (double q = 10; q < 1e4; q *= 9) {
    for (double x = 1e-7; x < 1; x *= 10) {
      for (int pid = -5; pid < 6; ++pid) {
        for (auto p : pdfs) {
          const double xf_g = p->xfxQ(pid, x, q);
          cout << pid << ", " << x << ", " << q << ": " << xf_g << endl;
        }
      }
    }
  }

  return 0;
}
