#include "bicubic.h"
#include "PDFGrid.h"

double ddx( const PDFGrid& Grid, const double* PID, const Idx_t XIdx, const Idx_t Q2Idx ) {    
    //Check for edge
    if (XIdx == 0) {
        //Use forward difference
        //return (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]) / (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]) * (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]);
        return (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]);
    }
    else if (XIdx == Grid.getXKnots().size()-1) {
	    //Use backward difference
        return (PID[Grid.index(XIdx, Q2Idx)] - PID[Grid.index(XIdx-1, Q2Idx)]);
    }
    else {
        //Use central difference
        double lddx = (PID[Grid.index(XIdx, Q2Idx)] - PID[Grid.index(XIdx-1, Q2Idx)]) / (Grid.getXKnots()[XIdx]-Grid.getXKnots()[XIdx-1]) * (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]);
        //double rddx = (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]) / (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]) * (Grid.getXKnots()[XIdx+1]-Grid.getXKnots()[XIdx]);
        double rddx = (PID[Grid.index(XIdx+1, Q2Idx)] - PID[Grid.index(XIdx, Q2Idx)]);
    
        return (0.5 * lddx + 0.5 * rddx);
    }
}

double interpolateCubic( const double T, const double YL, const double DL, const double YH, const double DH ) {
    //Pre calculate powers of T
    double t2 = T*T;
    double t3 = t2*T;
    
    //Calculate left point
    double p0 = (2*t3 - 3*t2 + 1)*YL;
    double m0 = (t3 - 2*t2 + T)*DL;
    
    //Calculate right point
    double p1 = (-2*t3 + 3*t2)*YH;
    double m1 = (t3 - t2)*DH;
    
    //Sum them
    return p0 + m0 + p1 + m1;
}

double BicubicInterpolator::interpolateQ2( const PDFGrid& Grid, const PID_t ID, const X_t X, const Q2_t Q2 ) const {
	//Index look up
	Idx_t xidx, q2idx;
	Grid.lookupClosestLow( X, Q2, xidx, q2idx );
	
	//Fractional parameter
	double tx = (X - Grid.getXKnots()[xidx]) / (Grid.getXKnots()[xidx+1] - Grid.getXKnots()[xidx]);
	double tq = (Q2 - Grid.getQ2Knots()[q2idx]) / (Grid.getQ2Knots()[q2idx+1] - Grid.getQ2Knots()[q2idx]);
	
	double* pid = Grid.getPIDData( ID );
	
	//Points in Q2
	double ql = interpolateCubic( tx, pid[Grid.index(xidx, q2idx)], ddx(Grid, pid, xidx, q2idx), pid[Grid.index(xidx+1, q2idx)], ddx(Grid, pid, xidx+1, q2idx) );
	double qh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+1)], ddx(Grid, pid, xidx, q2idx+1), pid[Grid.index(xidx+1, q2idx+1)], ddx(Grid, pid, xidx+1, q2idx+1) );
	
	//Derivatives in Q2
	double qdl, qdh;

	if (q2idx == 0) {
		//Forward difference for lower q
		qdl = (qh - ql);
		//qdl = (qh - ql) / (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]) * (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
		
		//Central difference for higher q
		double qhh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+1)], ddx(Grid, pid, xidx, q2idx+1), pid[Grid.index(xidx+1, q2idx+2)], ddx(Grid, pid, xidx+1, q2idx+2) );
		qdh = 0.5 * (qh - ql) + 0.5 * (qhh - qh)/(Grid.getQ2Knots()[q2idx+2]-Grid.getQ2Knots()[q2idx+1]) * (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
	}
	else if (q2idx == Grid.getQ2Knots().size()-1) {
	    //Central difference for lower q
	    double qll = interpolateCubic( tx, pid[Grid.index(xidx, q2idx-1)], ddx(Grid, pid, xidx, q2idx-1), pid[Grid.index(xidx+1, q2idx)], ddx(Grid, pid, xidx+1, q2idx) );
	    qdl = 0.5 * (qh - ql) + 0.5 * (ql - qll)/(Grid.getQ2Knots()[q2idx]-Grid.getQ2Knots()[q2idx-1]) * (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
        
        //Backward difference for higher q
        //qdh = (qh - ql);
        qdh = 0.0;
    }
    else {
        //Central difference for both q
        double qll = interpolateCubic( tx, pid[Grid.index(xidx, q2idx-1)], ddx(Grid, pid, xidx, q2idx-1), pid[Grid.index(xidx+1, q2idx)], ddx(Grid, pid, xidx+1, q2idx) );
        qdl = 0.5 * (qh - ql) + 0.5 * (ql - qll)/(Grid.getQ2Knots()[q2idx]-Grid.getQ2Knots()[q2idx-1]) * (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
                
        double qhh = interpolateCubic( tx, pid[Grid.index(xidx, q2idx+1)], ddx(Grid, pid, xidx, q2idx+1), pid[Grid.index(xidx+1, q2idx+2)], ddx(Grid, pid, xidx+1, q2idx+2) );
        qdh = 0.5 * (qh - ql) + 0.5 * (qhh - qh)/(Grid.getQ2Knots()[q2idx+2]-Grid.getQ2Knots()[q2idx+1]) * (Grid.getQ2Knots()[q2idx+1]-Grid.getQ2Knots()[q2idx]);
    }
    
    return interpolateCubic( tq, ql, qdl, qh, qdh );
}