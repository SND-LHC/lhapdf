//Include Guard
#pragma once

//LHAPDF Includes
#include "types.h"

namespace LHAPDF {
	//Forward Declarations
	class PDFGrid;
	
	//Extrapolator Declaration
	/// Extrapolator the general interface for extrapolating.
	class Extrapolator {
	public:
		/// Deconstruct Extrapolator.
		
		///
		virtual inline ~Extrapolator();
		
		/// Bind to PDFGrid.
		
		///
		/// \param Grid the PDFGrid to bind
		virtual inline void bind( const PDFGrid& );
		
		/// Unbind from PDFGrid.
		
		///
		/// \param Grid the PDFGrid to unbind
		virtual inline void unbind( const PDFGrid& );
	
		/// Interpolate single-point in Q.
		
		/// \param Grid the PDFGrid to extrapolate
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q the energy scale
		/// \return
		virtual inline double extrapolateQ( const PDFGrid&, const PID_t, const X_t, const Q_t ) const;
		
		/// Interpolator 
		
		/// \param Grid the PDFGrid to extrapolate
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return
		virtual double extrapolateQ2( const PDFGrid&, const PID_t, const X_t, const Q2_t ) const =0;
		
		///
		
		/// \param Grid
		/// \param ID the MC Parton
		/// \param PointNumber the number of points to extrapolate
		/// \param Xs 
		/// \param Qs
		/// \param XFXs
		virtual inline void extrapolateQ( const PDFGrid&, const PID_t, const uint32_t, const X_t*, const Q_t*, double* ) const;
		
		///
		
		/// \param Grid
		/// \param ID the MC Parton
		/// \param PointNumber the number of points to extrapolate
		/// \param Xs
		/// \param Q2s
		/// \param XFXs
		virtual inline void extrapolateQ2( const PDFGrid&, const PID_t, const uint32_t, const X_t*, const Q_t*, double* ) const;
	};
	
	//Extrapolator Defintions
	Extrapolator::~Extrapolator() {}
	void Extrapolator::bind( const PDFGrid& ) {}
	void Extrapolator::unbind( const PDFGrid& ) {}
	
	double Extrapolator::extrapolateQ( const PDFGrid& pdf, const PID_t id, const X_t x, const Q_t q ) const {
		return extrapolateQ2( pdf, id, x, q*q );
	}
	
	void Extrapolator::extrapolateQ( const PDFGrid& pdf, const PID_t id, const uint32_t numberofpoints, const X_t* xs, const Q_t* qs, double* samples ) const {
		//Applies individual extrapolation to all points in <std::vector<Point>&>
		for (uint32_t idx = 0; idx < numberofpoints; ++idx) {
			samples[idx] = extrapolateQ2( pdf, id, xs[idx], qs[idx]*qs[idx] );
		}
	}
	
	void Extrapolator::extrapolateQ2( const PDFGrid& pdf, const PID_t id, const uint32_t numberofpoints, const X_t* xs, const Q2_t* q2s, double* samples ) const {
		//Applies individual extrapolation to all points in <std::vector<Point>&>
		for (uint32_t idx = 0; idx < numberofpoints; ++idx) {
			samples[idx] = extrapolateQ2( pdf, id, xs[idx], q2s[idx] );
		}
	}
}