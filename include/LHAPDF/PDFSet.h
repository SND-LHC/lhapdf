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
      // std::cout << setinfopath.string() << std::endl;
      /// @todo Complain if the set info file has not been found? Is that mandatory?
      if (exists(setinfopath)) load(setinfopath);
      /// @todo Check that some mandatory metadata keys have been set? _check() function.
      /// @todo If not, try to guess some info or just exit?
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
      return metadata("SetDesc");
    }

    /// Get the type of PDF errors in this set (replica, symmhessian, asymmhessian, none)
    std::string errorType() const {
      return to_lower_copy(metadata("ErrorType"));
    }

    // /// Number of members in this set
    // int numMembers() const {
    //   return metadata_as<int>("NumMembers");
    // }

    // /// Number of members in this set
    size_t size() const {
      return metadata_as<unsigned int>("NumMembers");
    }

    //@}


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
    /// @note The supplied vector will be cleared before filling!
    ///
    /// @note As with the mkPDF method and factory function, the PDF pointers
    /// returned by this method are heap allocated and their memory management
    /// is now the responsibility of the caller.
    void mkPDFs(std::vector<PDF*>& pdfs) const {
      pdfs.clear();
      pdfs.reserve(size());
      for (size_t i = 0; i < size(); ++i) {
        pdfs.push_back( mkPDF(i) );
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

    //@}


    /// @todo Add AlphaS getter for set-level alphaS?


    /// @name Generic metadata cascading mechanism
    //@{

    /// Can this Info object return a value for the given key? (it may be defined non-locally)
    bool has_key(const std::string& key) const {
      return has_key_local(key) || getConfig().has_key(key);
    }

    /// Retrieve a metadata string by key name
    /// Rerite this cascading so that each stage only knows about the level above it
    const std::string& metadata(const std::string& key) const {
      if (has_key_local(key)) return metadata_local(key); //< value is defined locally
      return getConfig().metadata(key); //< fall back to the global config
    }

    //@}


  private:

    /// Name of this set
    std::string _setname;

  };


}
