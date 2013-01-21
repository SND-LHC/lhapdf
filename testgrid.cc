// Example program to test PDF grid format reading and interpolation

#include "LHAPDF/GridPDF.h"
#include <iostream>
using namespace LHAPDF;
using namespace std;

int main() {
  const GridPDF pdf("EXAMPLEPDF", 0);

  cout << "Verbosity = " << pdf.info().metadata("Verbosity") << endl;
  cout << "PdfDesc = " << pdf.info().metadata("PdfDesc") << endl;
  cout << "SetDesc = " << pdf.info().metadata("SetDesc") << endl;
  cout << "Flavors (str) = " << pdf.info().metadata("Flavors") << endl;
  vector<int> pids = pdf.info().metadata< vector<int> >("Flavors");
  cout << "Flavors (ints) = ";
  foreach (int f, pids) cout << f << " "; cout << endl;
  cout << "Flavors (vec<int>) = " << LHAPDF::to_str(pids) << endl;

  return 0;
}
