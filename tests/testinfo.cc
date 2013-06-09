// Example program for testing the info system

#include "LHAPDF/Info.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"
#include <iostream>
using namespace std;

int main() {
  const LHAPDF::Info cfg = LHAPDF::getConfig();
  cout << "UndefFlavorAction: " << cfg.metadata("UndefFlavorAction") << endl;
  cout << "Verbosity: " << cfg.metadata("Verbosity") << endl;

  const LHAPDF::PDFInfo info("CT10nlo", 2);

  if (info.has_key_local("PdfDesc")) {
    cout << 1.1 << endl;
  } else {
    cout << 1.2 << endl;
  }


  cout << 1.7 << endl;
  LHAPDF::getConfig();
  cout << 1.8 << endl;
  LHAPDF::getPDFSet("CT10nlo");
  cout << 1.9 << endl;


  cout << boolalpha << LHAPDF::getConfig().has_key("PdfDesc") << endl;
  cout << boolalpha << LHAPDF::getPDFSet("CT10nlo").has_key_local("PdfDesc") << endl;

  if (info.has_key_local("PdfType")) {
    cout << 2.1 << endl;
  } else {
    cout << 2.2 << endl;
  }

  if (info.has_key("PdfDesc")) {
    cout << 3.1 << endl;
    cout << "PdfDesc: " << info.metadata("PdfDesc") << endl;
  } else {
    cout << 3.2 << endl;
  }

  cout << "PdfType: " << info.metadata("PdfType") << endl;
  cout << "Verbosity: " << info.metadata("Verbosity") << endl;
  info.metadata_as<std::string>("PdfType");

  vector<int> pids = info.metadata_as< vector<int> >("Flavors");
  cout << "PIDs (1): "; foreach (int f, pids) { cout << f << " "; } cout << endl;
  cout << "PIDs (2): " << LHAPDF::to_str(pids) << endl;

  return 0;
}
