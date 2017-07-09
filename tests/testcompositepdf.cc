// Example program to test CompositePDFs

#include "LHAPDF/CompositePDF.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

  // Extract PDF ID strings from the command line
  if (argc < 2) cout << "Usage: testcompositepdf <PDF1> [<PDF2> ...]" << endl;
  vector<string> pdfstrs(&argv[1], &argv[argc]);
  if (pdfstrs.empty()) exit(1);

  // Create each requested PDF
  vector<LHAPDF::PDF*> pdfs;
  BOOST_FOREACH (const string& pstr, pdfstrs) {
    if (pstr.empty()) continue;
    LHAPDF::PDF* p = LHAPDF::mkPDF(pstr);
    pdfs.push_back(p);
  }

  // Create CompositePDF
  LHAPDF::CompositePDF cpdf(pdfs);
  cout << cpdf.xfxQ(21, 0.7, 10.0) << endl;
  cout << cpdf.xfxQ(21, 0.2, 126) << endl;
  /// @todo Ensure that PDFs have same flavour content? Return intersection of constituent PID sets?
  // BOOST_FOREACH (int pid, cpdf.flavors()) {
  //   cout << pid << " = " << cpdf.xfxQ(pid, 0.2, 124) << endl;
  // }

  return 0;
}
