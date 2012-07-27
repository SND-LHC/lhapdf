#include <iostream>
#include <math.h>

#include "LHAPDF/LHAPDF.h"

#include "PDFGrid.h"

int main(int argc, const char * argv[])
{
    //LHAPDF::initPDFSetByName( "MSTW2008lo90cl.LHgrid" );
    //LHAPDF::initPDF( 0 );
    
    /*srand( 9876 );
    
    std::string name("/Users/Martin/Development/LHAPDF/Tests/MSTW2008/MSTW2008lo90cl.LHm");
    PDFGrid* grid = PDFGrid::load( name );
    //PDFGrid* grid = PDFGrid::load( name, "bicubic" );
    
    //Interpolator* bic = createInterpolator( "bicubic" );
    //grid->setInterpolator( bic );
    
    uint32_t points = 10;
    uint32_t counter = 0;
    
    AxisPoints xs, qs;
    xs = grid->getXPoints();
    qs = grid->getQPoints();
    
    double xm, qm;
    double maxd = 0.0;
    double xfxmn, xfxmo;
    
    while( counter < points ) {
        double x, q;
        x = xs[0] + (xs[xs.size()-1] - xs[0])*((double)rand()/(double)RAND_MAX);
        q = qs[0] + (qs[qs.size()-1] - qs[0])*((double)rand()/(double)RAND_MAX);
        
        double xfxo = LHAPDF::xfx( x, q, 1 );
        double xfxn = grid->xfx( 1, x, q );
        
        if( fabs(xfxn - xfxo) > maxd )
        {
            xm = x;
            qm = q;
            xfxmn = xfxn;
            xfxmo = xfxo;
            
            maxd = fabs( xfxn - xfxo );
        }
        
        ++counter;
    }*/
    
    /*while( counter < points ) {
        double x, q;
        x = xs[0] + (xs[xs.size()-1] - xs[0])*((double)rand()/(double)RAND_MAX);
        q = qs[0] + (qs[qs.size()-1] - qs[0])*((double)rand()/(double)RAND_MAX);
        
        double xfxn = grid->xfx( 1, x, q );
        
        std::cout.width( 5 );
        std::cout.fill( ' ' );
        std::cout << counter;
        std::cout.width( 15 );
        std::cout.fill( ' ' );
        std::cout << x;
        std::cout.width( 15 );
        std::cout.fill( ' ' );
        std::cout << q;
        std::cout.width( 15 );
        std::cout.fill( ' ' );
        std::cout << xfxn << std::endl;
        
        ++counter;
    }*/
    
    return 0;
}
