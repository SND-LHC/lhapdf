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
    const std::string& metadata(const std::string& key) const;

    //@}


  private:

    std::string _setname;
    int _member;

  };


}
