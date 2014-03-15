/// Program to test LHAPDF6 automatic calculation of PDF uncertainties.
/// Written in March 2014 by G. Watt <Graeme.Watt(at)durham.ac.uk>.
/// Use formulae for PDF uncertainties and correlations in:
///   G. Watt, JHEP 1109 (2011) 069 [arXiv:1106.5788 [hep-ph]].

#include "LHAPDF/LHAPDF.h"
#include <boost/random.hpp>
using namespace std;

/// Simple test program to demonstrate the four PDFSet member functions.
///   set.uncertainty(values);
///   set.uncertainty(values, inputCL, median=false);
///   set.correlation(valuesA, valuesB);
///   set.randomValue(values, random, symmetrise=true);

int main(int argc, char* argv[]) {

  if (argc < 2) {
    cerr << "You must specify a PDF set:  ./testpdfunc setname" << endl;
    return 1;
  }

  const string setname = argv[1];
  const LHAPDF::PDFSet set(setname);
  const int nmem = set.size()-1;

  double x = 0.1; // momentum fraction
  double Q = 100.0; // factorisation scale in GeV

  /// Fill vectors xgAll and xuAll using all PDF members.
  /// Could replace xg, xu by cross section, acceptance etc.
  const vector<LHAPDF::PDF*> pdfs = set.mkPDFs();
  vector<double> xgAll, xuAll;
  for (int imem = 0; imem <= nmem; imem++) {
    xgAll.push_back(pdfs[imem]->xfxQ(21,x,Q)); // gluon distribution
    xuAll.push_back(pdfs[imem]->xfxQ(2,x,Q)); // up-quark distribution
  }

  /// Define formats for printing labels and numbers in output.
  string labformat = "%2s%10s%12s%12s%12s%12s\n";
  string numformat = "%12.4e%12.4e%12.4e%12.4e%12.4e\n";

  /// Calculate PDF uncertainty on gluon distribution.
  cout << "Gluon distribution at Q = " << Q << " GeV (normal uncertainties)" << endl;
  printf(labformat.c_str()," #","x","xg","error+","error-","error");
  vector<double> xgErr = set.uncertainty(xgAll);
  double xg0 = xgErr[0], xgp = xgErr[1], xgm = xgErr[2], xgs = xgErr[3];
  printf(numformat.c_str(),x,xg0,xgp,xgm,xgs);
  cout << endl;

  /// Calculate PDF uncertainty on up-quark distribution.
  cout << "Up-quark distribution at Q = " << Q << " GeV (normal uncertainties)" << endl;
  printf(labformat.c_str()," #","x","xu","error+","error-","error");
  vector<double> xuErr = set.uncertainty(xuAll);
  double xu0 = xuErr[0], xup = xuErr[1], xum = xuErr[2], xus = xuErr[3];
  printf(numformat.c_str(),x,xu0,xup,xum,xus);
  cout << endl;

  /// Calculate PDF correlation between gluon and up-quark.
  const double corr = set.correlation(xgAll,xuAll);
  cout << "Correlation between xg and xu = " << corr << endl;
  cout << endl;

  /// Calculate gluon PDF uncertainty scaled to 90% C.L.
  cout << "Gluon distribution at Q = " << Q << " GeV (scaled uncertainties)" << endl;
  printf(labformat.c_str()," #","x","xg","error+","error-","error");
  xgErr = set.uncertainty(xgAll, 90);
  xg0 = xgErr[0]; xgp = xgErr[1]; xgm = xgErr[2]; xgs = xgErr[3];
  printf(numformat.c_str(),x,xg0,xgp,xgm,xgs);
  cout << "Scaled PDF uncertainties to 90% C.L. using scale = " << xgErr[4] << endl;
  cout << endl;

  /// Calculate up-quark PDF uncertainty scaled to 1-sigma.
  /// Note: z-sigma = erf(z/sqrt(2)) = {0.68268949, 0.95449974, 0.99730020} for z = {1, 2, 3}.
  double sigma = 100*boost::math::erf(1/sqrt(2)); // 68.268949%
  cout << "Up-quark distribution at Q = " << Q << " GeV (scaled uncertainties)" << endl;
  printf(labformat.c_str()," #","x","xu","error+","error-","error");
  xuErr = set.uncertainty(xuAll, sigma);
  xu0 = xuErr[0]; xup = xuErr[1]; xum = xuErr[2]; xus = xuErr[3];
  printf(numformat.c_str(),x,xu0,xup,xum,xus);
  cout << "Scaled PDF uncertainties to 1-sigma using scale = " << xuErr[4] << endl;
  cout << endl;

  if (set.errorType() == "replicas") {

    /// Calculate gluon PDF as median and 90% C.L. of replica probability distribution.
    cout << "Gluon distribution at Q = " << Q << " GeV (median and 90% C.L.)" << endl;
    printf(labformat.c_str()," #","x","xg","error+","error-","error");
    xgErr = set.uncertainty(xgAll, 90, true);
    xg0 = xgErr[0]; xgp = xgErr[1]; xgm = xgErr[2]; xgs = xgErr[3];
    printf(numformat.c_str(),x,xg0,xgp,xgm,xgs);
    cout << endl;
    
    /// Calculate up-quark PDF as median and 68% C.L. of replica probability distribution.
    cout << "Up-quark distribution at Q = " << Q << " GeV (median and 68% C.L.)" << endl;
    printf(labformat.c_str()," #","x","xu","error+","error-","error");
    xuErr = set.uncertainty(xuAll, 68, true);
    xu0 = xuErr[0]; xup = xuErr[1]; xum = xuErr[2]; xus = xuErr[3];
    printf(numformat.c_str(),x,xu0,xup,xum,xus);
    cout << endl;
    
  } else if (set.errorType() == "hessian" || set.errorType() == "symmhessian") {

    /// Generate random values from Hessian best-fit and eigenvector values.
    /// See: G. Watt and R.S. Thorne, JHEP 1208 (2012) 052 [arXiv:1205.4024 [hep-ph]].
    /// Obtain Gaussian random numbers using Boost, but could also use C++11.

    int neigen = 0; // number of eigenvectors
    if (set.errorType() == "hessian") { 
      neigen = nmem/2;
    } else if (set.errorType() == "symmhessian") {
      neigen = nmem;
    }
    const unsigned seed = 1234;
    boost::mt19937 rng(seed);
    boost::normal_distribution<> nd;
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
    //C++11: default_random_engine generator(seed);
    //C++11: normal_distribution<double> distribution; // mean 0.0, s.d. = 1.0
    const int nrand = 5; // generate nrand random values
    for (int irand = 1; irand <= nrand; irand++) {
      /// Fill vector "random" with neigen Gaussian random numbers.
      vector<double> random;
      for (int ieigen=1; ieigen <= neigen; ieigen++) {
	double r = var_nor(); // using Boost
	//C++11: double r = distribution(generator); // using C++11
	random.push_back(r);
      }
      //const bool symmetrise = false; // average differs from best-fit
      const bool symmetrise = true; // default: average tends to best-fit
      double xgrand = set.randomValue(xgAll, random, symmetrise);
      /// Pass *same* random numbers to preserve correlations between xg and xu.
      double xurand = set.randomValue(xuAll, random, symmetrise);
      cout << "Random " << irand << ": xg = " << xgrand << ", xu = " << xurand << endl;
    }
    /// Random values generated in this way can subsequently be used for
    /// applications such as Bayesian reweighting or combining predictions
    /// from different groups (as an alternative to taking the envelope).
    /// See, for example, material at http://mstwpdf.hepforge.org/random/.

  }

  return 0;

}
