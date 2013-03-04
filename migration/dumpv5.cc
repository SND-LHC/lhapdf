#include <iostream>
#include <math.h>

#include <LHAPDF/LHAPDF.h>

int main( int argc, const char* argv[] ) {
    //Use CT10 LHgrid
    LHAPDF::initPDFSetByName( "CT10.LHgrid" );

    //Centeral value
    LHAPDF::initPDF( 0 );

    //Dump out points in (X,Q2)
    for( double x = 0.05; x <= 1.0; x += 0.05 ) {
        for( double q2 = 1e2; q2 < 1e5; q2 += 1e2 ) {
            std::cout << x << " " << q2 << " ";
            std::cout << LHAPDF::xfx( x, sqrt(q2), 2 ) << std::endl;
        }
    }
    
    return 0;
}
