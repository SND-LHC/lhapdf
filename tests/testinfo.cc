// Example program for testing the info system

#include "LHAPDF/Info.h"
#include <iostream>
using namespace std;

int main() {
  const LHAPDF::Info& cfg = LHAPDF::config();
  cout << "ImplicitFlavorAction: " << cfg.metadata("ImplicitFlavorAction") << endl;
  cout << "Verbosity: " << cfg.metadata("Verbosity") << endl;

  const LHAPDF::Info info("CT10nlo", 2);
  if (info.has_key("PdfDesc")) cout << "PdfDesc: " << info.metadata("PdfDesc") << endl;
  cout << "PdfType: " << info.metadata("PdfType") << endl;
  vector<int> pids = info.metadata< vector<int> >("Flavors");
  cout << "PIDs (1): "; foreach (int f, pids) { cout << f << " "; } cout << endl;
  cout << "PIDs (2): " << LHAPDF::to_str(pids) << endl;
  return 0;
}
