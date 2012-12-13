//LHAPDF Includes
#include "Nearest.h"
#include "PDFGrid.h"

namespace LHAPDF {
	//NearestExtrapolator Definitions
	double NearestExtrapolator::extrapolateQ2( const PDFGrid& grid, const PID_t id, const X_t x, const Q2_t q2 ) const {
		Idx_t xidx, q2idx;
		grid.lookupClosest( x, q2, xidx, q2idx );
		
		return grid.getPIDData(id)[grid.index(xidx, q2idx)];
	}
}
