#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <math.h>

#include <LHAPDF/LHAPDF.h>

#define MINLOGX -10
#define MINLOGQ2 1
#define MAXLOGQ2 8

#define DX 0.1
#define DQ2 0.1

int main( int argc, const char* argv[] ) {
    if( argc < 3 ) {
        std::cerr << "./dumpv5 name member." << std::endl;
        
        return 1;
    }

    //Use CT10 LHgrid
    const char* setname = argv[1];
    LHAPDF::initPDFSetByName( setname );

    //Centeral value
    const int member = atoi( argv[2] );
    LHAPDF::initPDF( member );

    //Dump out points in (X,Q2)
    for( int flavour = -6; flavour < 7; ++flavour ) {
        std::stringstream filename;
        filename << "flavour_v6_" << flavour << ".dat";

        std::ofstream output( filename.str().c_str() );
        for( double logX = MINLOGX; logX <= 0.0; logX += DX ) {
            for( double logQ2 = MINLOGQ2; logQ2 <= MAXLOGQ2; logQ2 += DQ2 ) {
                double x  = pow( 10, logX );
                double q2 = pow( 10, logQ2 );

                //std::cout << x << " " << q2 << " ";
                //std::cout << LHAPDF::xfx( x, sqrt(q2), flavour ) << std::endl;

                output << x << " " << q2 << " ";
                output << LHAPDF::xfx(x,sqrt(q2),flavour) << std::endl;
            }
        }

        output.close();
    }
    
    return 0;
}
