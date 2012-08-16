//Standard includes
#include <cstdio>
#include <iostream>
#include <math.h>
#include <cstdlib>

//LHAPDFv5 includes
#include "LHAPDF.h"

//LHAPDFv6 includes
#include "types.h"
#include "PDFGrid.h"
#include "ErrExtrapolator.h"

using namespace LHAPDF;

int main(int argc, const char * argv[]) {
    LHAPDF::initPDFSetByName( "MSTW2008lo90cl_nf3.LHgrid" );
    LHAPDF::initPDF( 0 );
    
    std::cout << "Loaded lhapdfv5..." << std::endl;

    std::string path("/Users/Martin/Desktop/lhapdfv6/tests/MSTW2008lo90cl_nf3/mbr_0.LHgm");
    //PDFGrid* grid = PDFGrid::load( path );
    //PDFGrid* grid = PDFGrid::load( path, "bicubic" );
    PDFSet* set = PDFSet::loadByName( "MSTW2008lo90cl_nf3" );
    std::cout << "Loaded set" << std::endl;
    PDFGrid* grid = dynamic_cast<PDFGrid*>( set->getMember( 0 ) );
    
    std::cout << "Loaded lhapdfv6..." << std::endl;
    
    AxisKnots xs, q2s;
    xs = grid->getXKnots();
    q2s = grid->getQ2Knots();
    
    AxisKnots::const_iterator c = xs.begin();
    for( ; c != xs.end(); ++c )
    	std::cout << *c << " ";
    std::cout << std::endl;
    
    ErrExtrapolator extra;
    grid->setExtrapolator( &extra );

	std::cout << "Starting..." << std::endl;

	//Accuracy
	for( unsigned int xidx = 0; xidx < xs.size()-1; ++xidx ) {
		for( unsigned int q2idx = 0; q2idx < q2s.size()-1; ++q2idx ) {
			//Interval middle
			double x = xs[xidx] * 0.5 + xs[xidx+1] * 0.5;
			double q2 = q2s[q2idx] * 0.5 + q2s[q2idx+1] * 0.5;
			
			double xfxv5 = LHAPDF::xfx( x, sqrt(q2), 1 );
			double xfxv6 = grid->xfxQ2( 1, x, q2 );
						
			double diff = fabs(xfxv6 - xfxv5);
			double p = fabs((xfxv5 - xfxv6) / xfxv5) * 100.0;
			
			if( p > 1.0 )
				printf( "@x %10.5e | @q2 %10.5e | @xfxv5 %15.5e | @xfxv6 %15.5e | %% %11f\n", x, q2, xfxv5, xfxv6, p );
		}
		
		std::cout << std::endl;
	}
	
	/*unsigned int xps = xs.size()*4, q2ps = q2s.size()*4;	
	for( unsigned int xp = 0; xp < xps-2; ++xp ) {
		for( unsigned int q2p = 0; q2p < q2ps-2; ++q2p ) {
			double x = 0.01 + xs[0] + (xs[xs.size()-1]-xs[0]) * (double)xp/(double)xps;
			double q2 = 0.01 + q2s[0] + (q2s[q2s.size()-1]-q2s[0]) * (double)q2p/(double)q2ps;
			
			double xfxv5 = LHAPDF::xfx( x, sqrt(q2), 2 );
			double xfxv6 = grid->xfxQ2( 2, x, q2 );
						
			double diff = fabs(xfxv6 - xfxv5);
			double p = fabs((xfxv5 - xfxv6) / xfxv5*100.0);
			
			if( p > 1.0 )
				printf( "@x %10.5f | @q2 %10.5f | %% %11f\n", x, q2, p );
		}
	}*/

	/*for( unsigned int q2idx = 0; q2idx < q2s.size()-1; ++q2idx ) {
		for( unsigned int xidx = 0; xidx < xs.size()-1; ++xidx ) {
			//Interval middle			
			double x = xs[xidx] * 0.5 + xs[xidx+1] * 0.5;
			double q2 = q2s[q2idx] * 0.5 + q2s[q2idx+1] * 0.5;
			
			double xfxv5 = LHAPDF::xfx( x, sqrt(q2), 2 );
			double xfxv6 = grid->xfxQ2( 2, x, q2 );
						
			double diff = fabs(xfxv6 - xfxv5);
			double p = fabs((xfxv5 - xfxv6) / xfxv5*100.0);
			
			if( p > 1.0 )
				printf( "@x %10e | @q2 %10e | %% %12f\n", x, q2, p );
		}
		
		std::cout << std::endl;
	}*/
	
	/*srand( 123456 );
	unsigned counter = 0;
	while (counter < 1000000) {
		double x = xs[0] + (xs[xs.size()-1]-xs[0]) * (double)rand()/(double)RAND_MAX;
		double q2 = q2s[0] + (q2s[q2s.size()-1]-q2s[0]) * (double)rand()/(double)RAND_MAX;
		
		double xfxv5 = LHAPDF::xfx( x, sqrt(q2), 2 );
		double xfxv6 = grid->xfxQ2( 2, x, q2 );
	
		++counter;
		
		//OUTPUT
		double diff = xfxv6 - xfxv5;
		double relerr = fabs(diff / xfxv5);
		double percent = relerr * 100.0;
		
		//if (percent > 1.0)
			printf( "@x %5.3e | @q2 %5.3e | %% %10f\n", x, q2, percent );
	}*/
	
    return 0;
}
