#include <iostream>
#include "PDF.h"
#include <cmath>
#include <iostream>

using namespace LHAPDF;
using namespace std;


/// Demo of an analytic/custom PDF class
struct PDFAnalytic: public PDF {

  double xfxQ2(PID_t id, double x, double q2) const {
    return 0.15 * sin(20.0*x) * sin(20.0*q2); // same for all PIDs
  }

  bool inRangeQ2(double x, const double q2) const {
    return true;
  }

  bool hasPID(double id) const {
    return true;
  }

};


int main( int argc, const char * argv[] ) {
  PDFAnalytic a;
  for (double x = 0; x < 1.0; x += 0.1) {
    for (double q2 = 0; q2 < 10; ++q2) {
      cout << x << " " << q2 << " " << a.xfxQ2( 0, x, q2 ) << endl;
    }
  }
  return 0;
}
