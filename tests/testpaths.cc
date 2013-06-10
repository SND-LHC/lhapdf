// Test program for path searching machinery

#include "LHAPDF/Paths.h"
#include <iostream>
using namespace std;

int main() {
  foreach (const LHAPDF::path& p, LHAPDF::paths()) {
    cout << p << endl;
  }
  cout << "@" << LHAPDF::findFile("lhapdf.conf") << "@" << endl;

  cout << "List of available PDFs:" << endl;
  foreach (const string& s, LHAPDF::availablePDFSets()) {
    cout << " " << s << endl;
  }

  return 0;
}
