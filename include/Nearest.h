//Inclusion Guards
#pragma once

//LHAPDF Includes
#include "types.h"
#include "Extrapolator.h"

//NearestExtrapolator Declaration
///Extrapolates nearest point of PDFGrid
class NearestExtrapolator: public Extrapolator
{
    /// Extrapolates using the closest point on the Grid.
    
    /// \param Grid the PDFGrid to extrapolate
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return the xfx value at (X,Q2)
    double extrapolateQ2( const PDFGrid&, const PID_t, const X_t, const Q2_t ) const;
};