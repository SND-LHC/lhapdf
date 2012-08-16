#include <iostream>
#include "PDFAnalytic.h"

using namespace LHAPDF;

int main( int argc, const char * argv[] ) {
    PDFAnalytic a;
    
    for( double x = 0; x < 1.0; x += 0.1 ) {
        for( double q2 = 0; q2 < 10; ++q2 ) {
            std::cout << x << " " << q2 << " " << a.xfxQ2( 0, x, q2 ) << std::endl;
        }
    }
    
    return 0;
}