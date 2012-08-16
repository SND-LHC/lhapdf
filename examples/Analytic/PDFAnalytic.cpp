#include "PDFAnalytic.h"
#include <cmath>

namespace LHAPDF {
	PDFAnalytic::PDFAnalytic() {
		//Does not construct anything...
	}
	
	PDFAnalytic::~PDFAnalytic() {
		//Has nothing to unload...
	}
	
	double PDFAnalytic::xfxQ2( const PID_t ID, const X_t X, const Q2_t Q2 ) const {
		//Same function for all pdfs	
		return 0.15*sin(20.0*X)*sin(20.0*Q2);
	}
	
	bool PDFAnalytic::inRangeQ2( const X_t X, const Q2_t Q2 ) const {
		return true;
	}
	
	bool PDFAnalytic::hasPID( const PID_t ID ) const {
		return true;
	}
}