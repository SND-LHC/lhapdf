#include "PDFAnalytic.h"
#include <cmath>

PDFAnalytic::PDFAnalytic() {
	//Generate discretized grid
	
	for( int x = 0; x < 10; ++x ) {
		for( int q = 0; q < 10; ++q ) {
			
		}
	}
}

PDFAnalytic::~PDFAnalytic() {
	
}

double PDFAnalytic::xfxQ2( const PID_t ID, const X_t X, const Q2_t Q2 ) const {
	
}

bool PDFAnalytic::inRangeQ2( const X_t X, const Q2_t Q2 ) const {
    return true;
}

bool PDFAnalytic::hasPID( const PID_t ID ) const {
    return true;
}