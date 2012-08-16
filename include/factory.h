#pragma once

//Standard Includes
#include <string>
#include <sstream>
#include <stdexcept>

///Include LHAPDF
#include "PDFSet.h"

///Include all Interpolator Types
#include "Bilinear.h"
#include "Bicubic.h"

///Include all Extrapolator Types
#include "Nearest.h"

namespace LHAPDF {
	///Interpolator Factory
	static Interpolator* createInterpolator( const std::string& name ) {
		//Insert interpolators here
		if (name == "bilinear")
			return new BilinearInterpolator();
		else if (name == "bicubic")
			return new BicubicInterpolator();
		else {
			std::stringstream error;
			error << "LHAPDFv6: Undeclared interpolator requested: " << name;
			
			throw std::runtime_error( error.str() );
		}
	}

	///Extrapolator Factory
	static Extrapolator* createExtrapolator( const std::string& name ) {
		//Insert extrapolators here
		if (name == "nearest")
			return new NearestExtrapolator();
		else {
			std::stringstream error;
			error << "LHAPDFv6: Undeclared extrapolator requested: " << name;
			
			throw std::runtime_error( error.str() );
		}
	}
}
