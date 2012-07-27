//Include Guard
#pragma once

//Standard Library Includes
#include <vector>
#include <exception>
#include <sstream>

//LHAPDF Includes
#include "types.h"

//PDF Declaration
/// PDF allows for any representation of PDFs.

/// The PDF interface declares the general form of all PDF types, such as Grid based or analytic.
class PDF {
public:
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
    
    //METADATA IMPLEMENTATION NEEDED!!!
    inline Member_t getMemberNumber() const;
    
    //inline access to meta data!!!
    //inline PDFSet* getPDFSet() comst;
    
    //DOES PDF NEED TO EXPOSE META DATA TO INTERPOLATORS??? YES, BMass, ZMass.... so we need to provide the meta data map, ALPAHS...
    
protected:
    //Metadata optimally map?? how can I store different data types in a strong type language
    uint32_t membernumber;
    //PDFSet* set; //This should store PID data for all members? set meta data or pdf meta data
};

//PDF Exception Declarations
/// UndefinedPIDException is thrown when an undefined PID_t is sampled from a PDF.
class UndefinedPIDException: public std::exception {
public:
    UndefinedPIDException( const PDF* PDF, const PID_t ID )
    : pdf(PDF), pid(ID) {}
    
    const char* what() const throw() {
        std::stringstream message;
        message << "Trying to access undefined PID_t " << pid << " in member " << pdf->getMemberNumber() << ".";
        return message.str().c_str();
    }
    
private:
    const PDF* pdf;
    const PID_t pid;
};

//PDF Definitions
bool PDF::inRangeQ( const X_t X, const Q_t Q ) const {
    return inRangeQ2( X, Q*Q );
}

double PDF::xfxQ( const PID_t ID, const X_t X, const Q_t Q ) const {
    return xfxQ2( ID, X, Q*Q );
}

void PDF::xfxQ( const PID_t ID, const uint32_t PointNumber, const X_t* Xs, const Q_t* Qs, double* XFXs ) const {
    //Check that PDF includes the PID_t requested
    if ( !hasPID( ID ) ) throw UndefinedPIDException( this, ID );
    
    //Iterate through all points given to get individual samples
    for( uint32_t idx = 0; idx < PointNumber; ++idx ) {
        //Call individual xfx per Point
        XFXs[idx] = xfxQ2( ID, Xs[idx], Qs[idx]*Qs[idx] );
    }
}

void PDF::xfxQ2( const PID_t ID, const uint32_t PointNumber, const X_t* Xs, const Q2_t* Q2s, double* XFXs ) const {
    //Check that PDF includes the PID_t requested
    if ( !hasPID( ID ) ) throw UndefinedPIDException( this, ID );
    
    //Iterate through all points given to get individual samples
    for( uint32_t idx = 0; idx < PointNumber; ++idx ) {
        //Call individual xfx per Point
        XFXs[idx] = xfxQ2( ID, Xs[idx], Q2s[idx] );
    }
}

Member_t PDF::getMemberNumber() const {
    return membernumber;
}