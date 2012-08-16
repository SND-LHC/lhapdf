//Include Guard
#pragma once

//LHAPDF Includes
#include "types.h"

namespace LHAPDF {
	//Forward Declarations
	class PDFGrid;
	
	//Interpolator Declaration
	/// Interpolator is the general interface for interpolating.
	class Interpolator
	{
	public:
		/// Deconstruct Interpolator.
		
		/// Default behaviour for Interpolator deconstruction is nothing.
		virtual inline ~Interpolator();
	
		/// Bind to PDFGrid. 
		
		/// Default behaviour for PDF binding is nothing.
		///
		/// \param Grid the PDFGrid to bind
		virtual inline void bind( const PDFGrid& );
		
		/// Unbind from PDFGrid.
		
		/// Default behaviour for PDF unbinding is nothing.
		///
		/// \param Grid the PDFGrid to unbind
		virtual inline void unbind( const PDFGrid& );
		
		/// Interpolate single-point in Q.
		
		/// \param Grid the PDFGrid to interpolate
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q the energy scale
		/// \return 
		virtual inline double interpolateQ( const PDFGrid&, const PID_t, const X_t, const Q_t ) const;
		
		/// Interpolate single-point in Q2.
		
		/// \param Grid the PDFGrid to interpolate
		/// \param ID the MC Parton ID
		/// \param X the momentum fraction
		/// \param Q2 the squared energy scale
		/// \return 
		virtual double interpolateQ2( const PDFGrid&, const PID_t, const X_t, const Q2_t ) const =0;
		
		/// Interpolate multi-point in Q.
		
		/// \param Grid the PDFGrid to interpolate
		/// \param ID the MC Parton ID
		/// \param PointNumber the number of points to interpolate
		/// \param Xs 
		/// \param Qs 
		/// \param XFXs 
		virtual inline void interpolateQ( const PDFGrid&, const PID_t, const uint32_t, const X_t*, const Q_t*, double* ) const;
		
		/// Interpolate multi-point in Q2.
		
		/// \param Grid the PDFGrid to interpolate
		/// \param ID the MC Parton ID
		/// \param PointNumber the number of points to interpolate
		/// \param Xs
		/// \param Q2s 
		/// \param XFXs the values of xfx at points (X,Q2) 
		virtual inline void interpolateQ2( const PDFGrid&, const PID_t, const uint32_t, const X_t*, const Q2_t*, double* ) const;
	};
	
	//Interpolator Definitions
	Interpolator::~Interpolator() {}
	void Interpolator::bind( const PDFGrid& pdf ) {}
	void Interpolator::unbind( const PDFGrid& pdf ) {}
	
	double Interpolator::interpolateQ( const PDFGrid& pdf, const PID_t id, const X_t x, const Q_t q ) const {
		return interpolateQ2( pdf, id, x, q*q );
	}
	
	void Interpolator::interpolateQ( const PDFGrid& pdf, const PID_t id, const uint32_t numberofpoints, const X_t* xs, const Q_t* qs, double* samples ) const {
		//Applies individual interpolation to all points in <std::vector<Point>&>
		for (uint32_t idx = 0; idx < numberofpoints; ++idx) {
			samples[idx] = interpolateQ2( pdf, id, xs[idx], qs[idx]*qs[idx] );
		}
	}
	
	void Interpolator::interpolateQ2( const PDFGrid& pdf, const PID_t id, const uint32_t numberofpoints, const X_t* xs, const Q2_t* q2s, double* samples ) const {
		//Applies individual interpolation to all points in <std::vector<Point>&>
		for (uint32_t idx = 0; idx < numberofpoints; ++idx) {
			samples[idx] = interpolateQ2( pdf, id, xs[idx], q2s[idx] );
		}
	}
}