// Test program for path searching machinery

#include "LHAPDF/Paths.h"
#include <iostream>

using namespace std;

int main() {
  foreach (const LHAPDF::path& p, LHAPDF::paths()) {
    cout << p << endl;
  }
  cout << "@" << LHAPDF::findFile("lhapdf.conf") << "@" << endl;
  return 0;
}
