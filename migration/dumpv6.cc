#include <iostream>

#include <LHAPDF/GridPDF.h>

int main( int argc, const char* argv[] ) {
    //Use CT10 PDFSet
    const LHAPDF::GridPDF pdf( "CT10", 0 );

    //Dump out points in (X,Q2)
    for( double x = 0.05; x <= 1.0; x += 0.05 ) {
        for( double q2 = 1e2; q2 < 1e5; q2 += 1e2 ) {
            std::cout << x << " " << q2 << " ";
            std::cout << pdf.xfxQ2( 2, x, q2 ) << std::endl;
        }
    }

    return 0;
}
