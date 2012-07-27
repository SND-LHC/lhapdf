#include "PDFGrid.h"

#include "Interpolator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

bool PDFGrid::hasInterpolator() const {
    if (interpolator != 0)
        return true;
    else 
        return false;
}

void PDFGrid::setInterpolator( Interpolator* I ) {
    if (interpolator != 0) {
        if (allocatedInterpolator)
            //Delete the allocated interpolator
            delete interpolator;
        else
            //Release association with PDFGrid
            interpolator->bind( *this );
    }
    
    //Set Interpolator
    interpolator = I;
    allocatedInterpolator = false;
    
    //Allow preprocessing of PDFGrid
    interpolator->bind( *this );
}

void PDFGrid::setDefaultInterpolator() {
    if (interpolator != 0) {
        if (allocatedInterpolator)
            //Delete the allocated interpolator
            delete interpolator;
        else
            //Release association with PDFGrid
            interpolator->unbind( *this );
    }
    
    //Use Interpolator factory
    interpolator = createInterpolator( "bilinear" );
    allocatedInterpolator = true;
}

bool PDFGrid::hasExtrapolator() const {
    if (extrapolator != 0)
        return true;
    else 
        return false;
}

void PDFGrid::setExtrapolator( Extrapolator* I ) {
    if (extrapolator != 0) {
        if (allocatedExtrapolator)
            delete extrapolator;
        else
            extrapolator->unbind( *this );
    }
    allocatedExtrapolator = false;
    
    extrapolator = I;
    extrapolator->bind( *this );
}

void PDFGrid::setDefaultExtrapolator() {
    if (extrapolator != 0) {
        if (allocatedExtrapolator)
            delete extrapolator;
        else
            extrapolator->unbind( *this );
    }

    //Use Extrapolator factory
    extrapolator = createExtrapolator( "nearest" );
    allocatedExtrapolator = true;
}

double PDFGrid::xfxQ2( const PID_t ID, const X_t X, const Q2_t Q2 ) const
{
    ///Error checks
    //Physical Ranges
    if (!inPhysicalRangeQ2( X, Q2 ))
        throw NonPhysicalRangeException( this, X, Q2 );
    //Undefined PIDs
    if (flavors.find( ID ) == flavors.end())
        throw UndefinedPIDException( this, ID );
    
    //OUTSIDE IS EXTRAPOLATION : INSIDE IS INTERPOLATION
    if (inRangeQ2( X, Q2 )) {
        //Valid Interpolator
        if( !hasInterpolator())
            throw UndefinedInterpolatorException( this );
        
        //Invoke interpolator to find value
        return interpolator->interpolateQ2( *this, ID, X, Q2 );
    }
    else {
        //Valid Extrapolator
        if (!hasExtrapolator())
            throw UndefinedExtrapolatorException( this );
        
        //Invoke extrapolator to find value
        return extrapolator->extrapolateQ2( *this, ID, X, Q2 );
    }
}

Idx_t bisectionSearch( const double X, const AxisKnots& Points ) {
    Idx_t low = 0, high = Points.size()-1;
    
    while ((high - low) > 1) {
        Idx_t middle = low + ( (high - low) >> 1 );
        
        if (X > Points[middle])
            low = middle;
        else
            high = middle;
    }
    
    return low;
}

void PDFGrid::lookupClosestLow( const X_t X, const Q2_t Q2, Idx_t& XIdx, Idx_t& Q2Idx ) const {    
    //DO ERROR CHECKS ON X AND Q?
    //TRUST INTERPOLATORS?
    
    //Bisection search to find closest low limit
    //POSSIBLE: Implement a hunting algorithm to have a guess
    XIdx = bisectionSearch( X, xknots );
    Q2Idx = bisectionSearch( Q2, q2knots );
}

Idx_t binarySearch( const double X, const AxisKnots& Points )
{
    Idx_t low = 0, high = Points.size()-1;
    
    while (true) {
        Idx_t middle = low + ((high - low) >> 2);
        
        if (middle == low || middle == high)
            break;
        
        if (X > Points[middle])
            low = middle;
        else
            high = middle;
    }
    
    return low;
}

void PDFGrid::lookupClosest( const X_t X, const Q2_t Q2, Idx_t& XIdx, Idx_t& QIdx ) const {
    //Binary Search Axes
    XIdx = binarySearch( X, xknots );
    QIdx = binarySearch( Q2, q2knots );
}

inline void tokenize( const std::string& Line, AxisKnots& Knots ) {
    size_t start = 0;
    size_t pos = 0;
    
    do {
        pos = Line.find_first_of( " ", start );
        
        Knots.push_back( atof( Line.substr( start, (pos - start) ).c_str() ) );
        
        start = pos+1;
    } while (pos != std::string::npos);
}

PDFGrid* PDFGrid::load( const std::string& Path ) {
    //File stream to member
    std::ifstream file( Path.c_str(), std::ifstream::in );
    
    //Initiate PDFGrid
    PDFGrid* grid = new PDFGrid();
    
    //Read header
    std::string line;
    std::stringstream header;
    
    while (true) {
        //Get next line
        getline( file, line );
        
        //Erase trailing space
        size_t found = line.find_last_not_of( " " );
        if (found != std::string::npos)
            line.erase(found+1);
        
        //Check for header end
        if (line == "---") {
            header << "---";
            break;
        }
        else 
            header << line << std::endl;
    }
    
    //Parse header
    //std::cout << header.str() << std::endl << std::endl;
    
    //Parse X line
    getline( file, line );
    tokenize( line, grid->xknots );
    
    //Parse Q line
    getline( file, line );
    tokenize( line, grid->q2knots );
    
    //Parse PID line
    getline( file, line );
    
    AxisKnots pids;
    tokenize( line, pids );
    
    //Allocate flavor data
    AxisKnots::const_iterator piditer = pids.begin();
    for( ; piditer != pids.end(); ++piditer ) {
        grid->flavors[*piditer] = new double[grid->xknots.size()*grid->q2knots.size()];
    }
    
    //Parse flavor data
    uint32_t lcount = 0;
    while (file.good()) {
        getline( file, line );
        
        size_t start = 0;
        size_t pos;
        AxisKnots::const_iterator piditer = pids.begin();
        for (; piditer != pids.end(); ++piditer) {
            pos = line.find_first_of( " ", start );
            
            grid->flavors[*piditer][lcount] = atof( line.substr( start, (pos - start) ).c_str() );
            
            start = pos+1;
        }
        
        ++lcount;
    }

    //Set defaults
    grid->setDefaultInterpolator();
    grid->setDefaultExtrapolator();
    
    return grid;
}

PDFGrid* PDFGrid::loadByName( const std::string& Name ) {
    //TODO: IMPLEMENT
    //Search for Set in LHAPATH... default directory
    //Load by path
    
    return 0;
}
