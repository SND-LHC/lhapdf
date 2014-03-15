// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Info.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Version.h"
#include "LHAPDF/Config.h"
/// @todo Keep this to the .cc files?
#include <boost/math/distributions/chi_squared.hpp>

namespace LHAPDF {


  // Forward declaration
  class PDF;


  /// Class for PDF set metadata and manipulation
  class PDFSet : public Info {
  public:

    /// @name Creation and deletion
    //@{

    /// Default constructor (for container compatibility)
    /// @todo Remove?
    PDFSet() { }

    /// Constructor from a set name
    /// @todo Remove?
    PDFSet(const std::string& setname) {
      _setname = setname;
      const string setinfopath = findpdfsetinfopath(setname);
      if (!file_exists(setinfopath))
        throw ReadError("Data file not found for PDF set '" + setname + "'");
      // /// Print out a banner if sufficient verbosity is enabled
      // const int verbosity = get_entry_as<int>("Verbosity", 1);
      // if (verbosity > 0) {
      //   std::cout << "Loading PDF set '" << setname << "'" << std::endl;
      //   print(std::cout, verbosity);
      // }
      // Load info file
      load(setinfopath);
      /// @todo Check that some mandatory metadata keys have been set: _check() function.
    }

    //@}


    /// @name PDF set metadata specialisations
    //@{

    /// @brief PDF set name
    ///
    /// @note _Not_ taken from the .info metadata file.
    std::string name() const {
      return _setname;
    }

    /// Description of the set
    std::string description() const {
      return get_entry("SetDesc");
    }

    /// First LHAPDF global index in this PDF set
    int lhapdfID() const {
      return get_entry_as<int>("SetIndex", -1);
    }

    /// Version of this PDF set's data files
    int dataversion() const {
      return get_entry_as<int>("DataVersion", -1);
    }

    /// Get the type of PDF errors in this set (replica, symmhessian, hessian, none)
    std::string errorType() const {
      return to_lower_copy(get_entry("ErrorType", "UNKNOWN"));
    }

    /// @brief Get the confidence level of the Hessian eigenvectors, in percent.
    ///
    /// If not defined, assume 1-sigma = erf(1/sqrt(2)) = 68.268949% by default.
    double errorConfLevel() const {
      return get_entry_as<double>("ErrorConfLevel", 100*boost::math::erf(1/sqrt(2)));
    }

    /// Number of members in this set
    // int numMembers() const {
    //   return get_entry_as<int>("NumMembers");
    // }
    size_t size() const {
      return get_entry_as<unsigned int>("NumMembers");
    }

    //@}


    /// Summary printout
    void print(std::ostream& os=std::cout, int verbosity=1) const {
      std::stringstream ss;
      if (verbosity > 0)
        ss << name() << ", version " << dataversion() << "; " << size() << " PDF members";
      if (verbosity > 1)
        ss << "\n" << description();
      os << ss.str() << std::endl;
    }


    /// @name Creating PDF members
    //@{


    /// Make the nth PDF member in this set, returned by pointer
    ///
    /// @note As with the mkPDF factory method, the PDF pointer returned by this
    /// method is heap allocated and its memory management is now the
    /// responsibility of the caller.
    PDF* mkPDF(int member) const {
      return LHAPDF::mkPDF(name(), member);
    }


    /// Make all the PDFs in this set, filling a supplied vector with PDF pointers
    ///
    /// This version may be preferred in many circumstances, since it can avoid
    /// the overhead of creating a new temporary vector.
    ///
    /// A vector of *smart* pointers can be used, for any smart pointer type which
    /// supports construction from a raw pointer argument, e.g. unique_ptr<PDF>(PDF*).
    ///
    /// @note The supplied vector will be cleared before filling!
    ///
    /// @note As with the mkPDF method and factory function, the PDF pointers
    /// returned by this method are heap allocated and their memory management
    /// is now the responsibility of the caller, either manually for raw pointers
    /// or automatically if smart pointers are used.
    ///
    /// @note Use an *appropriate* smart pointer, of course! This depends in detail on
    /// how you will use the PDF objects (do you want shared or unique pointers?), but
    /// they also need to be compatible with storage in STL containers. This is *not*
    /// the case with std::auto_ptr (which for this reason is replaced with
    /// std::unique_ptr in C++11).
    template <typename PTR>
    void mkPDFs(std::vector<PTR>& pdfs) const {
      const int v = verbosity();
      if (v > 0) {
        std::cout << "LHAPDF " << version() << " loading all " << size() << " PDFs in set " << name() << std::endl;
        this->print(std::cout, v);
      }
      pdfs.clear();
      pdfs.reserve(size());
      if (v < 2) setVerbosity(0); //< Disable every-member printout unless verbosity level is high
      for (size_t i = 0; i < size(); ++i) {
        pdfs.push_back( PTR(mkPDF(i)) );
      }
      setVerbosity(v);
    }

    /// Make all the PDFs in this set, returning as a vector of PDF pointers
    ///
    /// @note As with the mkPDF method and factory function, the PDF pointers
    /// returned by this method are heap allocated and their memory management
    /// is now the responsibility of the caller.
    std::vector<PDF*> mkPDFs() const {
      std::vector<PDF*> rtn;
      mkPDFs(rtn);
      return rtn;
    }

    /// @todo Use the following with default function template args if C++11 is being used
    // template <typename PTR=PDF*>
    template <typename PTR>
    std::vector<PTR> mkPDFs() const {
      std::vector<PTR> rtn;
      mkPDFs(rtn);
      return rtn;
    }

    //@}


    /// @todo Add AlphaS getter for set-level alphaS?


    /// @name Generic metadata cascading mechanism
    //@{

    /// Can this Info object return a value for the given key? (it may be defined non-locally)
    bool has_key(const std::string& key) const {
      return has_key_local(key) || getConfig().has_key(key);
    }

    /// Retrieve a metadata string by key name
    const std::string& get_entry(const std::string& key) const {
      if (has_key_local(key)) return get_entry_local(key); //< value is defined locally
      return getConfig().get_entry(key); //< fall back to the global config
    }

    /// Retrieve a metadata string by key name, with a fallback
    const std::string& get_entry(const std::string& key, const std::string& fallback) const {
      return Info::get_entry(key, fallback);
    }

    //@}


    /// @name PDF set uncertainty functions (written by G. Watt in March 2014)
    //@{

    /// Calculate central value and error with appropriate formulae using vector "values".
    /// If the PDF set is given in the form of replicas, the uncertainty is given by the
    /// standard deviation, and the central (average) value is not necessarily "values[0]"
    /// for quantities with a non-linear dependence on PDFs.  In the Hessian approach,
    /// the central value is the best-fit "values[0]" and the uncertainty is given by either
    /// the symmetric or asymmetric formula using eigenvector PDF sets.
    std::vector<double> uncertainty(const std::vector<double>& values) const {

      if (values.size() != size()) {
	throw UserError("Error in LHAPDF::PDFSet::uncertainty.  Input vector must contain values for all PDF members.");
      }

      /// The following five values will be returned in a vector<double>.
      double central = values[0], errplus = 0.0, errminus = 0.0, errsym = 0.0, scale = 1.0;
      const int nmem = size()-1;

      if (errorType() == "replicas") {

	/// Calculate the average and standard deviation
	/// using Eqs. (2.3) and (2.4) of arXiv:1106.5788v2.
	double av = 0.0, sd = 0.0;
	for (int imem = 1; imem <= nmem; imem++) {
	  av += values[imem];
	  sd += sqr(values[imem]);
	}
	av /= nmem; sd /= nmem;
	sd = nmem/(nmem-1.0)*(sd-av*av);
	if (sd > 0.0 && nmem > 1) sd = sqrt(sd);
	else sd = 0.0;
	central = av;
	errplus = errminus = errsym = sd;

      } else if (errorType() == "symmhessian") {

	for (int ieigen = 1; ieigen <= nmem; ieigen++) {
	  errsym += sqr(values[ieigen]-values[0]);
	}
	errsym = sqrt(errsym);
	errplus = errsym;
	errminus = errsym;

      } else if (errorType() == "hessian") {

	/// Calculate the asymmetric and symmetric Hessian uncertainties
	/// using Eqs. (2.1), (2.2) and (2.6) of arXiv:1106.5788v2.
	for (int ieigen = 1; ieigen <= nmem/2; ieigen++) {
	  errplus += sqr(max(max(values[2*ieigen-1]-values[0],values[2*ieigen]-values[0]),0.0));
	  errminus += sqr(max(max(values[0]-values[2*ieigen-1],values[0]-values[2*ieigen]),0.0));
	  errsym += sqr(values[2*ieigen-1]-values[2*ieigen]);
	}
	errsym = 0.5*sqrt(errsym);
	errplus = sqrt(errplus);
	errminus = sqrt(errminus);

      } else {

	throw MetadataError("\"ErrorType: " + errorType() + "\" not supported by LHAPDF::PDFSet::uncertainty.");

      }

      /// Return central value and asymmetric/symmetric PDF uncertainties as a vector<double>.
      std::vector<double> err;
      err.push_back(central);
      err.push_back(errplus);
      err.push_back(errminus);
      err.push_back(errsym);
      err.push_back(scale);
      return err;

    }


    /// Optional argument "inputCL" to rescale uncertainties to a particular confidence level.
    /// If the PDF set is given in the form of replicas, then optional argument "median" will
    /// instead calculate the median and confidence interval of the probability distribution.
    std::vector<double> uncertainty(const vector<double>& values, const double& inputCL, bool median=false) const {

      double scale = 1.0;
      std::vector<double> err = uncertainty(values);

      double reqCL = inputCL / 100.0; // convert from percentage
      double errCL = errorConfLevel() / 100.0; // convert from percentage

      /// Check that reqCl and errCL both lie between 0 and 1.
      if (reqCL < 0 || reqCL > 1 || errCL < 0 || errCL > 1) {
    	return err;
      }

      if (errorType() == "replicas" && median) {

    	/// Compute median and requested C.L. directly from probability distribution of replicas.

    	/// Sort "values" into increasing order, ignoring zeroth member (average over replicas).
	std::vector<double> sorted = values;
	std::sort(sorted.begin()+1,sorted.end());
    	const int nmem = size()-1;
    	/// Define central value to be median.
    	if (nmem%2) { // odd nmem => one middle value
    	  err[0] = sorted[nmem/2+1];
    	} else { // even nmem => average of two middle values
    	  err[0] = 0.5*(sorted[nmem/2] + sorted[nmem/2+1]);
    	}
    	/// Define uncertainties with a C.L. given by "reqCL".
    	int upper = round(0.5*(1+reqCL)*nmem); // round to nearest integer
    	int lower = 1+round(0.5*(1-reqCL)*nmem); // round to nearest integer
    	err[1] = sorted[upper] - err[0];
    	err[2] = err[0] - sorted[lower];
    	err[3] = 0.5*(err[1] + err[2]); // symmetrised

      } else {

    	/// Calculate the qth quantile of the chi-squared distribution with one degree of freedom.
    	/// Examples: quantile(dist, q) = {0.988946, 1, 2.70554, 3.84146, 4} for q = {0.68, 1-sigma, 0.90, 0.95, 2-sigma}.
    	boost::math::chi_squared dist(1);
    	double qerrCL = boost::math::quantile(dist,errCL);
    	double qreqCL = boost::math::quantile(dist,reqCL);
    	scale = sqrt(qreqCL/qerrCL);

    	for (int i = 1; i <= 3; i++) {
    	  err[i] *= scale;
    	}
    	err[4] = scale;
      }

      return err;

    }


    /// Calculate the PDF correlation using appropriate formulae given vectors "valuesA" and "valuesB".
    /// The correlation can vary between -1 and +1 where values close to {-1,0,+1} mean that the two
    /// quantities A and B are {anticorrelated,uncorrelated,correlated}, respectively.
    double correlation(const vector<double>& valuesA, const vector<double>& valuesB) const {

      if (valuesA.size() != size() || valuesB.size() != size()) {
	throw UserError("Error in LHAPDF::PDFSet::correlation.  Input vectors must contain values for all PDF members.");
      }

      double cor = 0.0;

      std::vector<double> errA = uncertainty(valuesA);
      std::vector<double> errB = uncertainty(valuesB);
      const int nmem = size()-1;

      if (errorType() == "replicas" && nmem > 1) {

    	/// Calculate the correlation using Eq. (2.7) of arXiv:1106.5788v2.

    	for (int imem = 1; imem <= nmem; imem++) {
    	  cor += valuesA[imem]*valuesB[imem];
    	}
    	cor = (cor/nmem - errA[0]*errB[0])/(errA[3]*errB[3]) * nmem/(nmem-1.0);

      } else if (errorType() == "symmhessian") {

    	for (int ieigen = 1; ieigen <= nmem; ieigen++) {
    	  cor += (valuesA[ieigen]-errA[0])*(valuesB[ieigen]-errB[0]);
    	}
    	cor /= errA[3]*errB[3];

      } else if (errorType() == "hessian") {

    	/// Calculate the correlation using Eq. (2.5) of arXiv:1106.5788v2.

    	for (int ieigen = 1; ieigen <= nmem/2; ieigen++) {
    	  cor += (valuesA[2*ieigen-1]-valuesA[2*ieigen])*(valuesB[2*ieigen-1]-valuesB[2*ieigen]);
    	}
    	cor /= 4.0*errA[3]*errB[3];

      }

      return cor;

    }


    /// Generate a random value from Hessian "values" and Gaussian random numbers.
    double randomValue(const vector<double>& values, const vector<double>& random, bool symmetrise=true) const {

      if (values.size() != size()) {
	throw UserError("Error in LHAPDF::PDFSet::randomValue.  Input vector must contain values for all PDF members.");
      }

      double frand = 0.0;

      double sigma = 100*boost::math::erf(1/sqrt(2));
      std::vector<double> err = uncertainty(values, sigma);
      double scale = err[4];
      int nmem = size()-1;

      /// Allocate number of eigenvectors based on ErrorType.
      unsigned neigen = 0;
      if (errorType() == "hessian") {
	neigen = nmem/2;
      } else if (errorType() == "symmhessian") {
	neigen = nmem;
      } else {
	throw UserError("Error in LHAPDF::PDFSet::randomValue.  This PDF set is not in the Hessian format.");
      }

      if (random.size() != neigen) {
	throw UserError("Error in LHAPDF::PDFSet::randomValue.  Input vector must contain random numbers for all eigenvectors.");
      }

      if (errorType() == "symmhessian") {

    	frand = values[0];
    	/// Loop over number of eigenvectors.
    	for (unsigned ieigen = 1; ieigen <= neigen; ieigen++) {
    	  double r = random[ieigen-1]; // Gaussian random number
    	  frand += r*abs(values[ieigen]-values[0])*scale;
    	}

      } else if (errorType() == "hessian") {

    	/// Use either Eq. (6.4) or Eq. (6.5) of arXiv:1205.4024v2.

    	frand = values[0];
    	/// Loop over number of eigenvectors.
    	for (unsigned ieigen = 1; ieigen <= neigen; ieigen++) {
    	  double r = random[ieigen-1]; // Gaussian random number
    	  if (symmetrise) {
    	    frand += 0.5*r*abs(values[2*ieigen-1]-values[2*ieigen])*scale;
    	  } else { // not symmetrised
    	    if (r<0.0) frand -= r*(values[2*ieigen]-values[0])*scale; // negative direction
    	    else frand += r*(values[2*ieigen-1]-values[0])*scale; // positive direction
    	  }
    	}

      }

      return frand;

    }

    //@}


  private:

    /// Name of this set
    std::string _setname;

  };


}
