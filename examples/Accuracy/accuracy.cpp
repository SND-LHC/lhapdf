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

    //PDFSet* set = PDFSet::loadByName( "MSTW2008lo90cl_nf3" );
    //PDFSet* set = PDFSet::load("/Users/Martin/Desktop/lhapdfv6/tests/constant/");
    PDFSet* set = PDFSet::load("/Users/Martin/Desktop/lhapdfv6/tests/curved/");
    std::cout << "Loaded set" << std::endl;
    PDFGrid* grid = dynamic_cast<PDFGrid*>( set->getMember( 0 ) );
    
    std::cout << "Loaded lhapdfv6..." << std::endl;
    
    AxisKnots xs, q2s;
    xs = grid->getXKnots();
    q2s = grid->getQ2Knots();
    
    std::cout << "SIZE: " << xs.size() << " " << q2s.size() << std::endl;
    
    //BilinearInterpolator bil;
    //grid->setInterpolator( &bil );
    
    ErrExtrapolator extra;
    grid->setExtrapolator( &extra );

	std::cout << "Starting..." << std::endl;
	
	const PIDdata data = grid->getPIDData( 1 );
	std::cout << "data" << std::endl;

	/*for( unsigned int idx = 0; idx < xs.size() * q2s.size(); ++idx ) {
		std::cout << data[idx] << std::endl;
	}*/

	//Accuracy
	uint32_t samples = 10;
	for( unsigned int xidx = 0; xidx < xs.size()-1; ++xidx ) {
		for( unsigned int q2idx = 0; q2idx < q2s.size()-1; ++q2idx ) {
	//for( unsigned int q2idx = 0; q2idx < q2s.size()-1; ++q2idx ) {
	//	for( unsigned int xidx = 0; xidx < xs.size()-1; ++xidx ) {
	
			double xi = (xs[xidx+1] - xs[xidx])/samples;
			double q2i = (q2s[q2idx+1] - q2s[q2idx])/samples;
	
			for( uint32_t x_sample = 0; x_sample < samples; ++x_sample ) {
				for( uint32_t q2_sample = 0; q2_sample < samples; ++q2_sample ) {
					//Interval middle
					double x = xs[xidx] + xi * x_sample;
					double q2 = q2s[q2idx] + q2i * q2_sample;
							
					//double x = xs[xidx];
					//double q2 = q2s[q2idx];
			
					//double xfxv5 = LHAPDF::xfx( x, sqrt(q2), 1 );
					//double xfxv5 = 1.0;
					double xfxv5 = x*x + q2*q2;
					
					double xfxv6 = grid->xfxQ2( 1, x, q2 );
						
					double diff = xfxv6 - xfxv5;
					double rel = diff / xfxv5;
					double p = rel * 100.0;
			
					//printf( "@x %+10.4e | @q2 %+10.4e | @file %+10.4e | @xfxv5 %+10.4e | @xfxv6 %+10.4e | D %+10.4e | R %+10.4e | %% %10.4f\n", x, q2, data[grid->index(xidx,q2idx)],xfxv5, xfxv6, diff, rel, p );
					printf( "@x %+10.4e | @q2 %+10.4e | @xfxv5 %+10.4e | @xfxv6 %+10.4e | D %+10.4e | R %+10.4e | %% %10.4f\n", x, q2,xfxv5, xfxv6, diff, rel, p );
					//printf( "@x %10.4e | @q2 %10.4e | @xfxv5 %10.4e | @xfxv6 %10.4e | Idx %10.4e\n", x, q2, xfxv5, xfxv6, data[grid->index(xidx,q2idx)] );
			
					//printf( "%10.4e %10.4e %10.4e %10.4e\n", x, q2, diff, rel );
				}
			}
		}
		
		std::cout << std::endl;
	}
	
	/*double x = 0.75;
	uint32_t samples = 100;
	for( uint32_t idx = 0; idx < q2s.size()-1; ++idx ) {
		double q2l = q2s[idx];
		double q2h = q2s[idx+1];
		double q2i = (q2h - q2l)/samples;
		
		for( uint32_t sample = 0; sample < samples; ++sample ) {
			double q2 = q2l + q2i * sample;
		
			double xfxv5 = fabs(LHAPDF::xfx( x, sqrt(q2), 1 ));
			double xfxv6 = fabs(grid->xfxQ2( 1, x, q2 ));
						
			double diff = xfxv6 - xfxv5;
			double rel = diff / xfxv5;
			double p = rel * 100.0;
			
			printf( "@x %+10.4e | @q2 %+10.4e | @xfxv5 %+10.4e | @xfxv6 %+10.4e | D %+10.4e | R %+10.4e | %% %10.4f\n", x, q2,xfxv5, xfxv6, diff, rel, p );
		}
	}*/
	
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
