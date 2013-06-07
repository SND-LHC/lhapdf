#include "LHAPDF/LHAPDF.h"

#include <string>
#include <vector>
#include <map>
using namespace std;

int main() {

  LHAPDF::PDFSet set("CT10nlo");
  vector<LHAPDF::PDF*> pdfs = set.mkPDFs();
  cout << set.size() << ", " << pdfs.size() << endl;

  foreach (LHAPDF::PDF* p, pdfs) {
    const double xf_g = p->xfxQ(21, 1e-3, 126.0);
    cout << xf_g << endl;
    delete p;
  }

  return 0;
}
