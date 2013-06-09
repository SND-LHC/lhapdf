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

    /// Default constructor (needed to store in a map)
    /// @todo Remove?
    PDFSet() { }

    /// Constructor from a set name
    /// @todo Remove?
    PDFSet(const std::string& setname) {
      _setname = setname;
      // const path setdirpath = findFile(setname);
      const path setinfopath = findpdfsetinfopath(setname);
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

    /// Number of members in this set
    int numMembers() const {
      return Info::metadata<int>("NumMembers");
    }

    /// A shorter, more STL-like alias for the number of members in this set
    size_t size() const {
      return numMembers();
    }

    //@}


    /// @name Creating PDF members
    //@{

    /// Make all the PDFs in this set, filling a supplied vector with PDF pointers
    ///
    /// This version may be preferred in many circumstances, since it can avoid
    /// the overhead of creating a new temporary vector.
    ///
    /// @note The supplied vector will be cleared before filling!
    ///
    /// @note As with the mkPDF factory method, the PDF pointers returned by
    /// this method are heap allocated and their memory management is now the
    /// responsibility of the caller.
    void mkPDFs(std::vector<PDF*>& pdfs) const {
      pdfs.clear();
      pdfs.reserve(numMembers());
      for (int i = 0; i < numMembers(); ++i) {
        pdfs.push_back( mkPDF(name(), i) );
      }
    }

    /// Make all the PDFs in this set, returning as a vector of PDF pointers
    ///
    /// @note As with the mkPDF factory method, the PDF pointers returned by
    /// this method are heap allocated and their memory management is now the
    /// responsibility of the caller.
    std::vector<PDF*> mkPDFs() const {
      std::vector<PDF*> rtn;
      mkPDFs(rtn);
      return rtn;
    }

    //@}


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
