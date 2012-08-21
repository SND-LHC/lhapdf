#include "Bicubic.h"
#include "PDFGrid.h"

#include <iostream>

namespace LHAPDF {
	// Provides d/dx at all grid locations
	double ddx( const PDFGrid& grid, const double* pid, const Idx_t xidx, const Idx_t q2idx ) {    
		//Check for edge
		if (xidx == 0) {
			//Use forward difference
			return (pid[grid.index(xidx+1, q2idx)] - pid[grid.index(xidx, q2idx)]) / (grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]);
		}
		else if (xidx == grid.getXKnots().size()-1) {
			//Use backward difference
			return (pid[grid.index(xidx, q2idx)] - pid[grid.index(xidx-1, q2idx)]) / (grid.getXKnots()[xidx]-grid.getXKnots()[xidx-1]);
		}
		else {
			//Use central difference
			double lddx = (pid[grid.index(xidx, q2idx)] - pid[grid.index(xidx-1, q2idx)]) / (grid.getXKnots()[xidx]-grid.getXKnots()[xidx-1]);
			double rddx = (pid[grid.index(xidx+1, q2idx)] - pid[grid.index(xidx, q2idx)]) / (grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]);
		
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
		double t = p0 + m0 + p1 + m1;
		return t;
	}
	
	double BicubicInterpolator::interpolateQ2( const PDFGrid& grid, const PID_t id, const X_t x, const Q2_t q2 ) const {
		//Index look up
		Idx_t xidx, q2idx;
		grid.lookupClosestLow( x, q2, xidx, q2idx );
		
		//Fractional parameters
		double tx = (x - grid.getXKnots()[xidx]) / (grid.getXKnots()[xidx+1] - grid.getXKnots()[xidx]);
		double tq = (q2 - grid.getQ2Knots()[q2idx]) / (grid.getQ2Knots()[q2idx+1] - grid.getQ2Knots()[q2idx]);
		
		double* pid = grid.getPIDData( id );
		
		//Points in Q2
		double vl = interpolateCubic( tx, pid[grid.index(xidx, q2idx)], 
										  ddx(grid, pid, xidx, q2idx)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]), 
										  pid[grid.index(xidx+1, q2idx)], 
										  ddx(grid, pid, xidx+1, q2idx)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
		double vh = interpolateCubic( tx, pid[grid.index(xidx, q2idx+1)], 
										  ddx(grid, pid, xidx, q2idx+1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]), 
										  pid[grid.index(xidx+1, q2idx+1)], 
										  ddx(grid, pid, xidx+1, q2idx+1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
		
		//Derivatives in Q2
		double vdl, vdh;
		
		if (q2idx == 0) {
			//Forward difference for lower q
			vdl = (vh - vl)/(grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]);
			
			//Central difference for higher q
			double vhh = interpolateCubic( tx, pid[grid.index(xidx, q2idx+2)], 
											   ddx(grid, pid, xidx, q2idx+2)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]), 
											   pid[grid.index(xidx+1, q2idx+2)], 
											   ddx(grid, pid, xidx+1, q2idx+2)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
			
			vdh = 0.5 * vdl + 0.5 * (vhh - vh)/(grid.getQ2Knots()[q2idx+2]-grid.getQ2Knots()[q2idx+1]);
		}
		else if (q2idx+1 == grid.getQ2Knots().size()-1) {
			//Backward difference for higher q
			vdh = (vh - vl)/(grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]);
		
			//Central difference for lower q
			double vll = interpolateCubic( tx, pid[grid.index(xidx, q2idx-1)], 
											   ddx(grid, pid, xidx, q2idx-1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]), 
											   pid[grid.index(xidx+1, q2idx-1)], 
											   ddx(grid, pid, xidx+1, q2idx-1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
											   
			vdl = 0.5 * vdh + 0.5 * (vl - vll)/(grid.getQ2Knots()[q2idx]-grid.getQ2Knots()[q2idx-1]);
		}
		else {
			//Central difference for both q
			double vll = interpolateCubic( tx, pid[grid.index(xidx, q2idx-1)], 
											   ddx(grid, pid, xidx, q2idx-1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]),
											   pid[grid.index(xidx+1, q2idx-1)], 
											   ddx(grid, pid, xidx+1, q2idx-1)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
											   
			vdl = 0.5 * (vh - vl)/(grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]) + 0.5 * (vl - vll)/(grid.getQ2Knots()[q2idx]-grid.getQ2Knots()[q2idx-1]);
					
			double vhh = interpolateCubic( tx, pid[grid.index(xidx, q2idx+2)], 
											   ddx(grid, pid, xidx, q2idx+2)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]), 
											   pid[grid.index(xidx+1, q2idx+2)], 
											   ddx(grid, pid, xidx+1, q2idx+2)*(grid.getXKnots()[xidx+1]-grid.getXKnots()[xidx]) );
											   
			vdh = 0.5 * (vh - vl)/(grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]) + 0.5 * (vhh - vh)/(grid.getQ2Knots()[q2idx+2]-grid.getQ2Knots()[q2idx+1]);
		}
		
		vdl *= (grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]);
		vdh *= (grid.getQ2Knots()[q2idx+1]-grid.getQ2Knots()[q2idx]);
		
		double f_f = interpolateCubic( tq, vl, vdl, vh, vdh );
		
		return f_f;
	}
}
