#include "AlphaS.h"
#include <cmath>
#include <vector>
#include <iostream>

static const double PI = 3.1415926535;

// Masses of quarks in GeV.  Used to calculate the number of quarks that are active at a given energy range Q2
static const double quarkMasses[] = {0.0017, 0.0041, 0.1, 1.29, 4.1, 172};

static double lambdaQCD = 0.4;

// Calculates the number of active quark flavours at energy scale Q2
int AlphaS::getNf( double Q2 ) {
    bool isNf = false;
    int nf = 0;
    while ( isNf == false ) {
	if ( quarkMasses[nf] <= Q2 && nf < 6) {
	    nf++;
	} else {
	    isNf = true;
	}
    }
    return nf;
}

// Calculates beta functions given the number of active flavours
std::vector<double> AlphaS::getBetas( int nf ) {
    std::vector<double> v;
    double beta0 = 11 - ((2 / 3.0) * nf);
    double beta1 = 51 - ((19 / 3.0) * nf);
    double beta2 = 2857 - ((5033 / 9.0) * nf) + ((325 / 27.0) * nf * nf);
    v.push_back(beta0);
    v.push_back(beta1);
    v.push_back(beta2);
    return v;
}

double AlphaS::meta( double Q2 ) {
    return 0.0;
}

double AlphaS::analytical( double Q2 ) {

    // Find the number of quarks with mass thresholds less than given energy scale Q2
    int nf = getNf( Q2 );

    // Calculate beta-functions
    std::vector<double> beta = getBetas( nf );

    // Pre-calculate powers of these 
    double beta02 = beta[0] * beta[0];
    double beta12 = beta[1] * beta[1];

    // Pre-calculate powers of ln(Q2/lambdaQCD) 
    double x = Q2 / lambdaQCD;
    double lnx = log(x);
    double lnlnx = log(lnx);
    double lnx2 = lnx * lnx;
    
    // Pre-calculate powers of expansion term y, where y is 1/ln(Q2/lambdaQCD)
    double y = 1 / lnx;
    double y2 = y * y;

    // Calculate expansion terms Ay[a_0*y^0 - a_1*y^1 + By^2[a_20 + a_21 - a_22]] where a_ij are series coefficients
    double A = 4 * PI / beta[0];
    double a_0 = 1;
    double a_1 = -2 * beta[1] * lnlnx / beta02;
    double B = 4 * beta12 / (beta02 * beta02);
    double a_20 = (lnlnx - 0.5) * (lnlnx - 0.5);
    double a_21 = beta[2] * beta[0] / (8 * beta12);
    double a_22 = 5 / 4.0;

    // Combine terms
    double alphaS = A * y * (a_0 - (a_1*y) + (B * y2 * (a_20 + a_21 - a_22)));
    
    return alphaS;

} 

// Calculates first order derivative, dy/dt, as it appears in differential equation 
double AlphaS::derivative( double t, double y, std::vector<double> beta ) {
    double d0 = beta[0]/(-2*PI);
    double d1 = beta[1]/(-4*PI*PI);
    double d2 = beta[2]/(-64*PI*PI*PI);
    double derivative = (d0*y*y) + (d1*y*y*y) + (d2*y*y*y*y);
    derivative /= t;
    return derivative;
}

// Solves differential equation in alphaS using an implementation of RK4;
double AlphaS::numerical( double alphaSmZ, double Q2 ) {

    // Step size
    double h = 0.001;

    // Fractional threshold to which we run in energy scale
    double accuracy = 0.1;

    // Mass of the Z-boson
    double m_z = 91.1876;

    // Check if we are going to run forward or backwards in energy scale towards target Q2 
    if ( Q2 < m_z ) {
	// Set h to step backwards if our target is less than where we are at
	h = (-1)*h;
    }

    // Number of active flavours used in beta function calculations
    double nf = getNf( Q2 );
    std::vector<double> beta = getBetas( nf );

    // Instantiate differential equation
    double y = alphaSmZ;
    double t = m_z;
    double f1, f2, f3, f4, k1, k2, k3, k4; 

    // Run in energy until we are within our defined threshold energy
    while ( (fabs(Q2 - t))/Q2 > accuracy ) {

	// Calculate RK4 increments 
	    
	// Increment based on the slope at the beginning of the interval using a simple Euler step
	f1 = derivative( t, y, beta );	
	k1 = h*f1;

	// Increments based on the slope at the midpoint of the interval
	f2 = derivative( t + (h/2.0), y + (k1/2.0), beta );
	k2 = h*f2;

	f3 = derivative( t + (h/2.0), y + (k2/2.0), beta );
	k3 = h*f3;

	// Increment based on the slope at the end of the interval
	f4 = derivative( t + h, y + k3, beta );
	k4 = h*f4;

	// Caclulate the weighted average of these increments and step forward in alphaS and energy scale
	y += (1/6.0)*(k1 + (2*k2) + (2*k3) + k4);
	t += h;

    }
    return y;

}

