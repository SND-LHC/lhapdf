#include "LHAPDF/GridPDF.h"
#include "LHAPDF/BicubicInterpolator.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>

int main( int argc, char** argv ) {
  GridPDF* grid = GridPDF::load("test.LHm");
  BicubicInterpolator bic;
  grid->setInterpolator( &bic );

  // srand(time(NULL));
  // for (int p = 0; p < 10000; ++p) {
  //   double x = (double)rand() / (double)RAND_MAX;
  //   double q2 = (double)rand() / (double)RAND_MAX;
  //   std::cout << x << " " << q2 << " " << grid->xfxQ2( 0, x, q2 ) << std::endl;
  // }

  for (int x = 0; x < 100; ++x) {
    double xd = x/100.0;
    for (int q2 = 0; q2 < 100; ++q2) {
      double q2d = q2/100.0;
      std::cout << xd << " " << q2d << " " << grid->xfxQ2( 0, xd, q2d ) << std::endl;
    }
  }

}
