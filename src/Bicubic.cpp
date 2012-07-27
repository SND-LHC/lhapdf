#include "bicubic.h"
#include "PDFGrid.h"

#include <iostream>

// Provides d/dx at all grid locations
double ddx( const PDFGrid& Grid, const double* PID, const Idx_t XIdx, const Idx_t Q2Idx ) {    
    //Check for edge
    if (XIdx == 0) {
        //Use forward difference
        return (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]) / (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]);
    }
    else if (XIdx == Grid.getXKnots().size()-1) {
    	//Use backward difference
		return (PID[Grid.index(XIdx, Q2Idx)] - PID[Grid.index(XIdx-1, Q2Idx)]) / (Grid.getXKnots()[XIdx]-Grid.getXKnots()[XIdx-1]);
    }
    else {
        //Use central difference
        double lddx = (PID[Grid.index(XIdx, Q2Idx)] - PID[Grid.index(XIdx-1, Q2Idx)]) / (Grid.getXKnots()[XIdx]-Grid.getXKnots()[XIdx-1]);
        double rddx = (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]) / (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]);
    
        return 0.5 * lddx + 0.5 * rddx;
    }
}

double interpolateCubic( const double T, const double VL, const double VDL, const double VH, const double VDH ) {
    //Pre calculate powers of T
    double t2 = T*T;
    double t3 = t2*T;
    
    //Calculate left point
    double p0 = (2*t3 - 3*t2 + 1)*VL;
    double m0 = (t3 - 2*t2 + T)*VDL;
    
    //Calculate right point
    double p1 = (-2*t3 + 3*t2)*VH;
    double m1 = (t3 - t2)*VDH;
    
    //Sum them
    return p0 + m0 + p1 + m1;
}

double BicubicInterpolator::interpolateQ2( const PDFGrid& Grid, const PID_t ID, const X_t X, const Q2_t Q2 ) const {
	//Index look up
	Idx_t xidx, q2idx;
	Grid.lookupClosestLow( X, Q2, xidx, q2idx );
	
	//Fractional parameters
	double tx = (X - Grid.getXKnots()[xidx]) / (Grid.getXKnots()[xidx+1] - Grid.getXKnots()[xidx]);
	double tq = (Q2 - Grid.getQ2Knots()[q2idx]) / (Grid.getQ2Knots()[q2idx+1] - Grid.getQ2Knots()[q2idx]);
	
	double* pid = Grid.getPIDData( ID );
	
	//Points in Q2
	double vl = interpolateCubic( tx, pid[Grid.index(xidx, q2idx)], 
									  ddx(Grid, pid, xidx, q2idx)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]), 
									  pid[Grid.index(xidx+1, q2idx)], 
									  ddx(Grid, pid, xidx+1, q2idx)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
	double vh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+1)], 
									  ddx(Grid, pid, xidx, q2idx+1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]), 
									  pid[Grid.index(xidx+1, q2idx+1)], 
									  ddx(Grid, pid, xidx+1, q2idx+1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
	
    //std::cout << " : ";
	//std::cout << "vl " << vl << " vh " << vh << " ";
		
	//Derivatives in Q2
	double vdl, vdh;
	
	if (q2idx == 0) {
		//Forward difference for lower q
		vdl = (vh - vl)/(Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
		
		//Central difference for higher q
		double vhh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+2)], 
										   ddx(Grid, pid, xidx, q2idx+2)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]), 
										   pid[Grid.index(xidx+1, q2idx+2)], 
										   ddx(Grid, pid, xidx+1, q2idx+2)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
		//std::cout << " vhh " << vhh;
		
		vdh = 0.5 * vdl + 0.5 * (vhh - vh)/(Grid.getQ2Knots()[q2idx+2]-Grid.getQ2Knots()[q2idx+1]);
	}
	else if (q2idx+1 == Grid.getQ2Knots().size()-1) {
        //Backward difference for higher q
        vdh = (vh - vl)/(Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
	
	    //Central difference for lower q
	    double vll = interpolateCubic( tx, pid[Grid.index(xidx, q2idx-1)], 
	    								   ddx(Grid, pid, xidx, q2idx-1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]), 
	    								   pid[Grid.index(xidx+1, q2idx-1)], 
	    								   ddx(Grid, pid, xidx+1, q2idx-1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
	    //std::cout << " vll " << vll; 
	    vdl = 0.5 * vdh + 0.5 * (vl - vll)/(Grid.getQ2Knots()[q2idx]-Grid.getQ2Knots()[q2idx-1]);
    }
    else {
        //Central difference for both q
        double vll = interpolateCubic( tx, pid[Grid.index(xidx, q2idx-1)], 
        								   ddx(Grid, pid, xidx, q2idx-1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]),
        								   pid[Grid.index(xidx+1, q2idx-1)], 
        								   ddx(Grid, pid, xidx+1, q2idx-1)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
        //std::cout << " vll " << vll;
        vdl = 0.5 * (vh - vl)/(Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]) + 0.5 * (vl - vll)/(Grid.getQ2Knots()[q2idx]-Grid.getQ2Knots()[q2idx-1]);
                
        double vhh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+2)], 
        								   ddx(Grid, pid, xidx, q2idx+2)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]), 
        								   pid[Grid.index(xidx+1, q2idx+2)], 
        								   ddx(Grid, pid, xidx+1, q2idx+2)*(Grid.getXKnots()[xidx+1]-Grid.getXKnots()[xidx]) );
    	//std::cout << " vhh " << vhh;
        vdh = 0.5 * (vh - vl)/(Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]) + 0.5 * (vhh - vh)/(Grid.getQ2Knots()[q2idx+2]-Grid.getQ2Knots()[q2idx+1]);
    }
    
    //std::cout << " vdl " << vdl << " vdh " << vdh;
    //std::cout << " : ";
    
    vdl *= (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
    vdh *= (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
    
    return interpolateCubic( tq, vl, vdl, vh, vdh );
}