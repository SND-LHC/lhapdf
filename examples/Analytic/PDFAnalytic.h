#pragma once

#include "PDF.h"

class PDFAnalytic: public PDF {
public:
	PDFAnalytic();
	~PDFAnalytic();
	
    double xfxQ2( const PID_t, const X_t, const Q2_t ) const;
    
    bool inRangeQ2( const X_t, const Q2_t ) const;
    
    bool hasPID( const PID_t ) const;
};