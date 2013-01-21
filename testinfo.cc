// Example program for testing the info system

#include "LHAPDF/Info.h"
#include <iostream>
using namespace std;

int main() {
  const LHAPDF::Info& cfg = LHAPDF::Info::config();
  cout << cfg.metadata("ImplicitFlavorAction") << endl;
  cout << cfg.metadata("Verbosity") << endl;

  const LHAPDF::Info info("EXAMPLEPDF/EXAMPLEPDF_0000.lha");
  cout << info.metadata("PdfDesc") << endl;
  cout << info.metadata("PdfType") << endl;
  vector<int> pids = info.metadata< vector<int> >("Flavors");
  foreach (int f, pids) cout << f << endl;
  cout << LHAPDF::to_str(pids) << endl;
  return 0;
}
