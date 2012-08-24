//Include Guard
#pragma once

//Standard Library Includes
#include <vector>
#include <sstream>
#include <stdexcept>
#include <map>

//LHAPDF Includes
#include "types.h"

namespace LHAPDF {
	class PDFSet;
	
	//PDF Declaration
	/// PDF allows for any representation of PDFs.

	/// The PDF interface declares the general form of all PDF types, such as Grid based or analytic.
	class PDF {
	public:
		
		virtual inline ~PDF();
		
		/// Single-point xfx value at (X,Q) in PDF
		
		/// xfxQ will square the given Q and return the value from xfxQ2. All grids are defined to be in
		/// Q2 since this is the general case.
		///
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q the energy scale
		/// \return the value of xfx at (X,Q)
		virtual inline double xfxQ( const PID_t, const X_t, const Q_t ) const;
		
		/// Single-point xfx value at (X,Q2) in PDF
		
		/// xfxQ2 will return the interpolated or extrapolated value at (x,Q2) from the PDF.
		///
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q2 the energy scale squared
		/// \return the value of xfx at (X,Q2)
		virtual double xfxQ2( const PID_t, const X_t, const Q2_t ) const =0;
		
		/// Multi-point xfx values at given (x,Q) in PDF
		
		/// The multi-point xfx allows for parallel arrays of (X,Q) coordinates to be passed and will
		/// return an array as xfx values.
		///
		/// \param ID the MC Parton ID
		/// \param PointNumber the number of points
		/// \param Xs the momentum fraction of each point
		/// \param Qs the energy scale of each point
		/// \param XFXs the values of xfx at (X,Q)
		virtual inline void xfxQ( const PID_t, const uint32_t, const X_t*, const Q_t*, double* ) const;
		
		/// Multi-point xfx values at given (X,Q2) in PDF
	
		/// The multi-point xfx allows for parallel arrays of (X,Q) coordinates to be passed and will
		/// return an array as xfx values.
		///
		/// \param ID the MC Parton ID
		/// \param PointNumber the number of points
		/// \param Xs the momentum fraction of each point
		/// \param Qs the squared energy scale of each point
		/// \param XFXs the values of xfx at (X,Q2)
		virtual inline void xfxQ2( const PID_t, const uint32_t, const X_t*, const Q2_t*, double* ) const;
		
		/// Range check for (X,Q)
	
		/// inRangeQ will return true when given (X,Q) are in the range of this PDF.
		/// It actually squares given Q and returns value from inRangeQ2.
		///
		/// \param X the momentum fraction
		/// \param Q the energy scale
		/// \return 
		virtual inline bool inRangeQ( const X_t, const Q_t ) const;
		
		/// Range check for (X,Q2)
	
		/// inRangeQ2 will return true when given (X,Q2) are in the range of this PDF.
		///
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return 
		virtual bool inRangeQ2( const X_t, const Q2_t ) const =0;
		
		/// Parton component check
	
		/// Checks whether the ID is a valid parton in the PDF.
		///
		/// \param ID the MC Parton ID
		/// \return 
		virtual bool hasPID( const PID_t ) const =0;
		
		/// Gives reference to containing PDFSet
		inline const PDFSet& getPDFSet() const;
		
		/// Gives access to any meta data in the member
		inline std::string getMetadata( const std::string& ) const;
		
		/// Gives Member Name
		inline std::string getMemberName() const;
		
		/// Gives Member ID
		inline Member_t getMemberID() const;

		/// Gives access to the type it declares
		inline std::string getType() const;
		
	protected:
		const PDFSet* set;
		
		std::map<std::string, std::string> metadict;
	};
	
	//PDF Definitions
	PDF::~PDF() {}
	
	bool PDF::inRangeQ( const X_t x, const Q_t q ) const {
		return inRangeQ2( x, q*q );
	}
	
	double PDF::xfxQ( const PID_t id, const X_t x, const Q_t q ) const {
		return xfxQ2( id, x, q*q );
	}
	
	void PDF::xfxQ( const PID_t id, const uint32_t PointNumber, const X_t* xs, const Q_t* qs, double* samples ) const {
		//Check that PDF includes the PID_t requested
		if ( !hasPID( id ) ) {
			std::stringstream error;
			error << "Undefined PID_t requested: " << id;
		
			throw std::runtime_error( error.str() );
		}
		
		//Iterate through all points given to get individual samples
		for( uint32_t idx = 0; idx < PointNumber; ++idx ) {
			//Call individual xfx per Point
			samples[idx] = xfxQ2( id, xs[idx], qs[idx]*qs[idx] );
		}
	}
	
	void PDF::xfxQ2( const PID_t id, const uint32_t pointcount, const X_t* xs, const Q2_t* q2s, double* samples ) const {
		//Check that PDF includes the PID_t requested
		if ( !hasPID( id ) ) {
			std::stringstream error;
			error << "Undefined PID_t requested: " << id;
		
			throw std::runtime_error( error.str() );
		}		
		//Iterate through all points given to get individual samples
		for( uint32_t idx = 0; idx < pointcount; ++idx ) {
			//Call individual xfx per Point
			samples[idx] = xfxQ2( id, xs[idx], q2s[idx] );
		}
	}
	
	std::string PDF::getMemberName() const {
		return metadict.find("name")->second;
	}

	Member_t PDF::getMemberID() const {
		std::stringstream idss;
		idss << metadict.find("id")->second;
		
		Member_t id;
		idss >> id;
		
		return id;
	}
	
	std::string PDF::getType() const {
		return metadict.find("type")->second;
	}
	
	std::string PDF::getMetadata( const std::string& key ) const {
		std::map<std::string, std::string>::const_iterator data = metadict.find( key );
		
		if( data != metadict.end() )
			return data->second;
		else {
			//Metadata not found
			std::stringstream error;
			error << "Metadata for key: " << key << " not found.";
			
			throw std::runtime_error( error.str() );
		}
	}
}