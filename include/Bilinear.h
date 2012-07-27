//Include Guard
#pragma once

//LHAPDF Includes
#include "Interpolator.h"

/// The Bilinear interpolation implementation of Interpolator
class BilinearInterpolator: public Interpolator
{
public:
    /// 
    
    /// \param Grid the PDFGrid
    /// \param ID the MC Parton ID
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    double interpolateQ2( const PDFGrid&, const PID_t, const X_t, const Q2_t ) const;
};