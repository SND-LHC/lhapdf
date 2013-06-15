// Example program for testing the info system

#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"
#include <iostream>
using namespace std;

int main() {

  LHAPDF::Info& cfg = LHAPDF::getConfig();
  cout << "UndefFlavorAction: " << cfg.metadata("UndefFlavorAction") << endl;
  cout << "Verbosity: " << cfg.metadata("Verbosity") << endl;
  cfg.setMetadata("Verbosity", 5);
  const LHAPDF::Info& cfg2 = LHAPDF::getConfig();
  cout << "New Verbosity from second Config: " << cfg2.metadata("Verbosity") << endl;

  const LHAPDF::PDFSet set("CT10nlo");
  cout << "SetDesc: " << set.metadata("SetDesc") << endl;
  cout << "Verbosity from set: " << set.metadata("Verbosity") << endl;

  const LHAPDF::PDFInfo info("CT10nlo", 2);
  if (info.has_key("PdfDesc")) cout << "PdfDesc: " << info.metadata("PdfDesc") << endl;
  cout << "PdfType: " << info.metadata("PdfType") << endl;
  cout << "Verbosity from PDF: " << info.metadata("Verbosity") << endl;
  vector<int> pids = info.metadata_as< vector<int> >("Flavors");
  cout << "PIDs (1): "; foreach (int f, pids) { cout << f << " "; } cout << endl;
  cout << "PIDs (2): " << LHAPDF::to_str(pids) << endl;

  return 0;
}
