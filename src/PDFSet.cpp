#include "PDFSet.h"
#include <math.h>

static const double PI = 3.1415926535;
static const double quarkMasses[] = {0.0017, 0.0041, 0.1, 1.29, 4.1, 172};

// TODO: maybe allow for general lambdaQCD, set to 0.4 for now
static double lambdaQCD = 0.4;

PDFSet* PDFSet::load( const std::string& ) {
    //Load meta data for set
    
    return 0;
}

PDFSet* PDFSet::loadByName( const std::string& Name ) {
    //Search default path for name
    
    //Use load by path to load file
    
    return load( Name );
}

double PDFSet::alphaS( double Q2 ) {

    // Find the number of quarks with mass thresholds less than given energy scale Q2
    bool isNf = false;
    int nf = 0;
    while ( isNf == false ) {
	if ( quarkMasses[nf] <= Q2 ) {
	    nf++;
	} else {
	    isNf = true;
	}
    }

    // Calculate beta-functions
    double beta0 = 11 - ((2 / 3.0) * nf);
    double beta1 = 51 - ((19 / 3.0) * nf);
    double beta2 = 2857 - ((5033 / 9.0) * nf) + ((325 / 27.0) * nf * nf);

    // Pre-calculate powers of these 
    double beta02 = beta0 * beta0;
    double beta12 = beta1 * beta1;

    // Pre-calculate powers of ln(Q2/lambdaQCD) 
    double x = Q2 / lambdaQCD;
    double lnx = log(x);
    double lnlnx = log(lnx);
    double lnx2 = lnx * lnx;
    
    // Pre-calculate powers of expansion term y, where y is 1/ln(Q2/lambdaQCD)
    double y = 1 / lnx;
    double y2 = y * y;

    // Calculate expansion terms Ay[a_0*y^0 - a_1*y^1 + By^2[a_20 + a_21 - a_22]] where a_ij are series coefficients
    double A = 4 * PI / beta0;
    double a_0 = 1;
    double a_1 = -2 * beta1 * lnlnx / beta02;
    double B = 4 * beta12 / (beta02 * beta02);
    double a_20 = (lnlnx - 0.5) * (lnlnx - 0.5);
    double a_21 = beta2 * beta0 / (8 * beta12);
    double a_22 = 5 / 4.0;

    // Combine terms
    double alphaS = A * y * (a_0 - (a_1*y) + (B * y2 * (a_20 + a_21 - a_22)));
    
    return alphaS;

}  
