// Example program to test PDF grid format reading and interpolation

#include "LHAPDF/PDFGrid.h"
#include <iostream>
using namespace LHAPDF;
using namespace std;

int main() {
  const PDFGrid pdf("EXAMPLEPDF", 0);

  // cout << pdf.info().metadata("Foo") << endl;
  cout << pdf.info().metadata("Verbosity") << endl;
  cout << pdf.info().metadata("PdfDesc") << endl;
  cout << pdf.info().metadata("SetDesc") << endl;
  cout << pdf.info().metadata("Flavors") << endl;
  vector<int> pids = pdf.info().metadata< vector<int> >("Flavors");
  foreach (int f, pids) cout << f << " "; cout << endl;
  cout << LHAPDF::to_str(pids) << endl;

  return 0;
}
