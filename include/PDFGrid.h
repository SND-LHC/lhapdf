//Include Guards
#pragma once

//Standard Library Includes
#include <vector>
#include <map>
#include <stdlib.h>
#include <exception>
#include <fstream>

//LHAPDF Includes
#include "types.h"
#include "factory.h"
#include "PDF.h"
#include "PDFSet.h"

//YAML Includes
#include "yaml-cpp/yaml.h"

namespace LHAPDF {
	//Forward Declarations
	class Interpolator;
	class Extrapolator;
	
	//PDFGrid Declaration
	/// Grid implementation of PDF interface.
	class PDFGrid: public PDF {
	public:
		///
		inline PDFGrid( const PDFSet* );
	
		/// Deletes allocated data.
		inline ~PDFGrid();
		
		/// xfx value at (X,Q2).
		
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return 
		double xfxQ2( const PID_t, const X_t, const Q2_t ) const;
		
		/// Checks whether (X,Q2) is a valid physical point.
		
		/// This will return false in the following cases:
		///  1) on points below 0.0 in X, since a negative momentum fraction is not valid.
		///  2) on points above 1.0 in X, since X is a fraction of momentum.
		///  3) on points below 0.0 in Q2, TODO
		///
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return 
		inline bool inPhysicalRangeQ2( const X_t, const Q2_t ) const;
		
		/// Checks whether (X,Q2) is in PDF data.
		
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return 
		inline bool inRangeQ2( const X_t, const Q2_t ) const;
		
		/// Checks whether Parton exists in this PDF.
	
		/// \param PID Parton identification
		/// \return a boolean of whether the Parton is supported.
		inline bool hasPID( const PID_t ) const;
		
		/// Looks up the lowest knot of the rectangle in which (X,Q2) is located.
	
		/// This function is useful for interpolators to be able to look up which patch the (X,Q2) is
		/// in and then interpolating between the returned index and the index+1.
		///
		/// The index returned will never be on the high edge of the grid knots.
		///
		/// \param X Momentum fraction coordinate.
		/// \param Q2 Energy scale coordinate.
		/// \param XIdx Index in X axis.
		/// \param Q2Idx Index in Q2 axis.
		void lookupClosestLow( const X_t, const Q2_t, Idx_t&, Idx_t& ) const;
		
		/// Looks up the closest grid knot.
	
		/// \param X Momentum fraction
		/// \param Q2 Energy scale
		/// \param XIdx Index in X axis.
		/// \param Q2Idx Index in Q2 axis.
		void lookupClosest( const X_t, const Q2_t, Idx_t&, Idx_t& ) const;
		
		/// Checks whether the (XIdx,Q2Idx) is valid in this PDF.
		bool isValidIndex( const Idx_t, const Idx_t ) const;
		
		/// Transforms a (XIdx, Q2Idx) into a 1D index to access flavor images at knots in this PDF.
		inline Idx_t index( const Idx_t, const Idx_t ) const;
		
		/// Allows access to the X grid knots.
		inline const AxisKnots& getXKnots() const;
		
		/// Allows access to the Q2 grid knots.
		inline const AxisKnots& getQ2Knots() const;
		
		/// Allows access to any flavor image.
		inline const PIDdata getPIDData( const PID_t ) const;
		
		/// Sets Interpolator to be used for interpolating between grid knots.
		void setInterpolator( Interpolator* );
		
		/// Sets the default interpolator.
		void setDefaultInterpolator();
		
		/// Sets Extrapolator to be used for extrapolating outside grid knots.
		void setExtrapolator( Extrapolator* );
		
		/// Sets the default extrapolator.
		void setDefaultExtrapolator();
		
		/// Check of valid Interpolator.
		bool hasInterpolator() const;
		
		/// Check of valid Extrapolator.
		bool hasExtrapolator() const;
	
		/// Loads the given member by path
		
		/// \param Path the file path to the ".LHm" file
		/// \return 
		//static PDFGrid* load( const std::string&, const PDFSet& set );
		static PDFGrid* load( PDFGrid*, const YAML::Node&, std::ifstream& );
		
	private:		
		/// Momentum fraction knots
		AxisKnots xknots;
		
		/// Squared energy scale knots
		AxisKnots q2knots;
		
		/// Flavor Map
		PIDmap flavors;
		
		/// Associated Interpolator
		Interpolator* interpolator;
		/// Flag for self allocated Interpolator
		bool allocatedInterpolator;
		
		/// Associated Extrapolator
		Extrapolator* extrapolator;
		/// Flag for self allocated Extrapolator
		bool allocatedExtrapolator;
	};
		
	//PDFGrid Definitions
	PDFGrid::PDFGrid( const PDFSet* setp )
	: interpolator(0), allocatedInterpolator(false),
	  extrapolator(0), allocatedExtrapolator(false) {
		set = setp;  	
	}
	
	PDFGrid::~PDFGrid() {
		//Check if the class created the interpolator
		if (interpolator != 0 && allocatedInterpolator)
			delete interpolator;
		if (extrapolator != 0 && allocatedExtrapolator)
			delete extrapolator;
	}
	
	bool PDFGrid::inRangeQ2( const X_t x, const Q2_t q2 ) const {
		//Check x is in range
		if (x < xknots[0] || (xknots[xknots.size()-1] < x) )
			return false;
		//Check q2 is in range
		if (q2 < q2knots[0] || (q2knots[q2knots.size()-1] < q2) )
			return false;
		
		return true;
	}
	
	bool PDFGrid::inPhysicalRangeQ2( const X_t x, const Q2_t q2 ) const {
		//Check X is physical
		if( x < 0.0 || x > 1.0 )
			return false;
		//Check Q2 is physical
		else if( q2 < 0.0 )
			return false;
		else
			return true;
	}
	
	bool PDFGrid::hasPID( const PID_t id ) const {
		//Check flavor exists
		if (flavors.find( id ) == flavors.end())
			return false;
		else
			return true;
	}
	
	size_t PDFGrid::index( const Idx_t xidx, const Idx_t q2idx ) const {
		//X is along a scanline, Q2 is a scanline
		return xidx + q2idx * xknots.size();
	}
	
	const AxisKnots& PDFGrid::getXKnots() const {
		//Return X Knots
		return xknots;
	}
	
	const AxisKnots& PDFGrid::getQ2Knots() const {
		//Return Q2 Knots
		return q2knots;
	}
	
	const PIDdata PDFGrid::getPIDData( const PID_t id ) const {
		//Check PID exists
		if (!hasPID( id )) {
			std::stringstream error;
			error << "Undefined PID_t requested: " << id;
			
			throw std::runtime_error( error.str() );
		}
		
		return flavors.find( id )->second;
	}
}