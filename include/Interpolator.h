//Include Guard
#pragma once

//LHAPDF Includes
#include "types.h"

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
void Interpolator::bind( const PDFGrid& PDF ) {}
void Interpolator::unbind( const PDFGrid& PDF ) {}

double Interpolator::interpolateQ( const PDFGrid& PDF, const PID_t ID, const X_t X, const Q_t Q ) const {
    return interpolateQ2( PDF, ID, X, Q );
}

void Interpolator::interpolateQ( const PDFGrid& PDF, const PID_t ID, const uint32_t NumberOfPoints, const X_t* Xs, const Q_t* Qs, double* XFXs ) const {
    //Applies individual interpolation to all points in <std::vector<Point>&>
    for (uint32_t idx = 0; idx < NumberOfPoints; ++idx) {
        XFXs[idx] = interpolateQ2( PDF, ID, Xs[idx], Qs[idx]*Qs[idx] );
    }
}

void Interpolator::interpolateQ2( const PDFGrid& PDF, const PID_t ID, const uint32_t NumberOfPoints, const X_t* Xs, const Q2_t* Q2s, double* XFXs ) const {
    //Applies individual interpolation to all points in <std::vector<Point>&>
    for (uint32_t idx = 0; idx < NumberOfPoints; ++idx) {
        XFXs[idx] = interpolateQ2( PDF, ID, Xs[idx], Q2s[idx] );
    }
}