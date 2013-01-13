// Test program for path searching machinery

#include "LHAPDF/Paths.h"
#include <iostream>

using namespace std;

int main() {
  cout << 1 << endl;
  foreach (const boost::filesystem::path& p, LHAPDF::paths()) {
    cout << p << endl;
  }
  cout << 2 << endl;
  cout << "@" << LHAPDF::findFile("lhapdf.conf") << "@" << endl;
  return 0;
}
