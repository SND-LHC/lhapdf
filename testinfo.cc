// Example program for testing the info system

#include "LHAPDF/Info.h"
#include <iostream>
using namespace std;

int main() {
  const LHAPDF::Info& cfg = LHAPDF::Info::config();
  cout << cfg.metadata("Foo") << endl;
  cout << cfg.metadata("Flavors") << endl;
  vector<int> pids = cfg.metadata< vector<int> >("Flavors");
  foreach (int f, pids) cout << f << endl;
  cout << LHAPDF::to_str(pids) << endl;
  return 0;
}
