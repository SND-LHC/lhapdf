#include "LHAPDF/PDF.h"
#include <cmath>
#include <iostream>

using namespace std;


/// Demo of an analytic/custom PDF class
struct PDFAnalytic : public LHAPDF::PDF {

  double _xfxQ2(int id, double x, double q2) const {
    if (abs(id) > 5) return 0;
    return 0.15 * sin(20.0*x) * sin(20.0*q2); // same for all light flavors
  }

  bool inRangeX(double x) const { return true; }
  bool inRangeQ2(double q2) const { return true; }

  vector<int> flavors(double id) const {
    vector<int> rtn;
    for (int id = 1; id < 6; ++id) {
      rtn.push_back(id);
      rtn.push_back(-id);
    }
    sort(rtn.begin(), rtn.end());
    return rtn;
  }

};


int main( int argc, const char * argv[] ) {
  PDFAnalytic a;
  for (double x = 0; x < 1.0; x += 0.1) {
    for (double q2 = 0; q2 < 10; ++q2) {
      cout << x << " " << q2 << " " << a.xfxQ2(0, x, q2) << endl;
    }
  }
  return 0;
}
