//LHAPDF Includes
#include "nearest.h"
#include "PDFGrid.h"

//NearestExtrapolator Definitions
double NearestExtrapolator::extrapolateQ2( const PDFGrid& Grid, const PID_t ID, const X_t X, const Q2_t Q2 ) const {
    Idx_t xidx, q2idx;
    Grid.lookupClosest( X, Q2, xidx, q2idx );
    
    return Grid.getPIDData(ID)[Grid.index(xidx, q2idx)];
}