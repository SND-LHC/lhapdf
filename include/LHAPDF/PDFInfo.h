// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Info.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/PDFIndex.h"

namespace LHAPDF {


  /// Metadata class for PDF members
  class PDFInfo : public Info {
  public:

    /// @name Creation and deletion
    //@{

    /// Default constructor (for container compatibility)
    /// @todo Remove?
    PDFInfo() { }

    /// Constructor from a set name and member ID.
    PDFInfo(const std::string& setname, int member);

    /// Constructor from an LHAPDF ID code.
    PDFInfo(int lhaid);

    //@}


    /// @name Metadata accessors
    //@{

    /// Can this Info object return a value for the given key? (it may be defined non-locally)
    bool has_key(const std::string& key) const;

    /// Retrieve a metadata string by key name
    const std::string& get_entry(const std::string& key) const;

    /// Retrieve a metadata string by key name, with a fallback
    const std::string& get_entry(const std::string& key, const std::string& fallback) const {
      return Info::get_entry(key, fallback);
    }

    //@}


  private:

    /// Name of the set in which this PDF is contained (for PDFSet lookup)
    std::string _setname;

    /// Member ID in PDF set
    /// @note Not currently used, but could be useful if a memberID method is exposed.
    int _member;

  };


}
