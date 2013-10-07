// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Info.h"
#include "LHAPDF/Factories.h"

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
      const path setinfopath = findpdfsetinfopath(setname);
      if (!exists(setinfopath))
        throw ReadError("Data file not found for PDF set '" + setname + "'");
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

    /// Version of this PDF set's data files
    int dataversion() const {
      return get_entry_as<int>("DataVersion", -1);
    }

    /// Get the type of PDF errors in this set (replica, symmhessian, asymmhessian, none)
    std::string errorType() const {
      return to_lower_copy(get_entry("ErrorType", "UNKNOWN"));
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
    void print(std::ostream& os=std::cout) const {
      os << name() << ", version " << dataversion() << "\n"
         << size() << " PDF members\n"
         << description() << std::endl;
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
      pdfs.clear();
      pdfs.reserve(size());
      for (size_t i = 0; i < size(); ++i) {
        pdfs.push_back( PTR(mkPDF(i)) );
      }
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


  private:

    /// Name of this set
    std::string _setname;

  };


}
