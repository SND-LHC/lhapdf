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

namespace LHAPDF {


  // Forward declaration
  class PDF;


  /// Structure for storage of uncertainty info calculated over a PDF error set
  /// @todo Exact role of scale and default value?
  struct PDFUncertainty {
    PDFUncertainty(double cent=0, double eplus=0, double eminus=0, double esymm=0, double scalefactor=1)
      : central(cent), errplus(eplus), errminus(eminus), errsymm(esymm), scale(scalefactor)
    {    }

    double central, errplus, errminus, errsymm, scale;

    /// Cast operator to a vector<double>? Hmm...
    operator std::vector<double> () const {
      vector<double> rtn;
      rtn.resize(5);
      rtn.push_back(this->central);
      rtn.push_back(this->errplus);
      rtn.push_back(this->errminus);
      rtn.push_back(this->errsymm);
      rtn.push_back(this->scale);
      return rtn;
    }
  };


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
    PDFSet(const std::string& setname);

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
    double errorConfLevel() const;

    /// Number of members in this set
    // int numMembers() const {
    //   return get_entry_as<int>("NumMembers");
    // }
    size_t size() const {
      return get_entry_as<unsigned int>("NumMembers");
    }

    //@}


    /// Summary printout
    void print(std::ostream& os=std::cout, int verbosity=1) const;


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
    //
    /// @todo Needs to be implemented in the header since the arg type is templated.
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


    /// @name PDF set uncertainty functions
    //@{

    /// @brief Calculate central value and error with appropriate formulae using vector "values".
    ///
    /// If the PDF set is given in the form of replicas, the uncertainty is given by the
    /// standard deviation, and the central (average) value is not necessarily "values[0]"
    /// for quantities with a non-linear dependence on PDFs.  In the Hessian approach,
    /// the central value is the best-fit "values[0]" and the uncertainty is given by either
    /// the symmetric or asymmetric formula using eigenvector PDF sets.
    PDFUncertainty uncertainty(const std::vector<double>& values) const;

    /// @brief Calculate central value and error using vector @c values and optional args for stat treatment
    ///
    /// @todo Unify under this signature
    ///
    /// Optional argument @c inputCL to rescale uncertainties to a particular confidence level.
    /// If the PDF set is given in the form of replicas, then optional argument @c median will
    /// instead calculate the median and confidence interval of the probability distribution.
    PDFUncertainty uncertainty(const std::vector<double>& values, double inputCL, bool median=false) const;

    /// Calculate PDF uncertainties (as above), with with efficient no-copy return to the @c rtn argument.
    void uncertainty(PDFUncertainty& rtn, const std::vector<double>& values, double inputCL, bool median=false) const {
      rtn = uncertainty(values, inputCL, median);
    }

    /// @brief Calculate the PDF correlation between @c valuesA and @c valuesB using appropriate formulae for this set.
    ///
    /// The correlation can vary between -1 and +1 where values close to {-1,0,+1} mean that the two
    /// quantities A and B are {anticorrelated,uncorrelated,correlated}, respectively.
    double correlation(const vector<double>& valuesA, const vector<double>& valuesB) const;

    /// @brief Generate a random value from Hessian @c values and Gaussian random numbers.
    ///
    /// @todo Behaviour if this is not a Hessian set?
    double randomValue(const vector<double>& values, const vector<double>& random, bool symmetrise=true) const;

    //@}


  private:

    /// Name of this set
    std::string _setname;

  };


}
