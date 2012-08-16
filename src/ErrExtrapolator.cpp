//Standard Includes
#include <iostream>

//LHAPDF Includes
#include "ErrExtrapolator.h"
#include "PDFGrid.h"

namespace LHAPDF {
	//ErrExtrapoltor Definitions
	double ErrExtrapolator::extrapolateQ2( const PDFGrid& grid, const PID_t id, const X_t x, const Q2_t q2 ) const {
		std::cout << "Extrapolation not allowed." << std::endl;
		
		throw -1;
	}
}