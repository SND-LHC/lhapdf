#include <iostream>
#include "PDFAnalytic.h"

int main( int argc, const char * argv[] ) {
    PDFAnalytic a;
    
    for( int x = 0; x < 10; ++x ) {
        for( int q = 0; q < 10; ++q ) {
            double xd = x / 10.0;
            double qd = q / 10.0;
            
            std::cout << xd << " " << qd << " " << a.xfxQ( 0, xd, qd ) << std::endl;
        }
    }
    
    return 0;
}