#include <iostream>
#include <sstream>
#include <fstream>

#include <LHAPDF/GridPDF.h>

#define MINLOGX -8
#define MINLOGQ2 1
#define MAXLOGQ2 8

#define DX 0.1
#define DQ2 0.1

int main( int argc, const char* argv[] ) {
    if( argc < 3 ) {
        std::cerr << "./dumpv6 setname setmember" << std::endl;

        return 1;
    }

    //Use CT10 PDFSet
    const char* setname = argv[1];
    const int setmember = atoi( argv[2] );
    const LHAPDF::GridPDF pdf( setname, setmember );

    //Dump out points in (X,Q2)
    int flavours[] ={-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
    for( int idx = 0; idx < 11; ++idx ) {
        int flavour = flavours[idx];

        std::stringstream filename;
        filename << "flavour_v6_" << flavour << ".dat";

        std::ofstream output( filename.str().c_str() );
        for( double logX = MINLOGX; logX <= 0.0; logX += DX ) {
            for( double logQ2 = MINLOGQ2; logQ2 <= MAXLOGQ2; logQ2 += DQ2 ) {
                double x  = pow( 10, logX );
                double q2 = pow( 10, logQ2 );

                output << x << " " << q2 << " ";
                output << pdf.xfxQ2(flavour, x, q2) << std::endl;
            }
        }

        output.close();
    }

    return 0;
}
