// Example program for testing the info system

#include "LHAPDF/Info.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/Factories.h"
#include <iostream>
using namespace std;

int main() {
  const LHAPDF::Info cfg = LHAPDF::getConfig();
  cout << "UndefFlavorAction: " << cfg.metadata("UndefFlavorAction") << endl;
  cout << "Verbosity: " << cfg.metadata("Verbosity") << endl;

  const LHAPDF::PDFInfo info("CT10nlo", 2);
  if (info.has_key("PdfDesc")) cout << "PdfDesc: " << info.metadata("PdfDesc") << endl;
  cout << "PdfType: " << info.metadata("PdfType") << endl;
  cout << "Verbosity: " << info.metadata("Verbosity") << endl;
  info.metadata_as<std::string>("PdfType");

  vector<int> pids = info.metadata_as< vector<int> >("Flavors");
  cout << "PIDs (1): "; foreach (int f, pids) { cout << f << " "; } cout << endl;
  cout << "PIDs (2): " << LHAPDF::to_str(pids) << endl;

  return 0;
}
