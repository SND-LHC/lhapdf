// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_CompositePDF_H
#define LHAPDF_CompositePDF_H

#include "LHAPDF/PDF.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/Extrapolator.h"
#include "LHAPDF/KnotArray.h"

namespace LHAPDF {


  /// @brief A PDF defined by multiplicative combination of several constituent PDFs
  ///
  /// @todo Add division? And other operators? Becomes a can of worms quite fast...
  class CompositePDF : public PDF {
  public:

    /// @name Creation and deletion
    //@{

    /// Default (empty) constructor
    CompositePDF() {  }

    /// @brief Constructor from a list of PDF pointers
    ///
    /// Pointers must be heap-allocated (new'd), and ownership will be taken by the new object.
    template <typename PDFPTR>
    CompositePDF(const std::vector<PDFPTR>& pdfs) {
      setConstituentPDFs(pdfs);
    }

    /// Constructor from a list of set names and member IDs
    CompositePDF(const std::vector< std::pair<std::string, int> >& setnames_members);

    /// Constructor from a list of LHAPDF IDs
    CompositePDF(const std::vector<int>& lhaids);

    /// Virtual destructor to allow inheritance
    virtual ~CompositePDF() {
      reset();
    }

    //@}


    /// @name Handling the constituent PDFs
    //@{

    /// Get the list of constituent PDFs (const version)
    const std::vector<PDF*> constituentPDFs() const {
      return _pdfs;
    }

    /// Get the list of constituent PDFs (non-const version)
    std::vector<PDF*> constituentPDFs() {
      return _pdfs;
    }

    /// Set the list of constituent PDFs
    template <typename PDFPTR>
    void setConstituentPDFs(const std::vector<PDFPTR>& pdfs) {
      reset();
      addConstituentPDFs(pdfs);
    }

    /// Append a PDF to the list of constituents
    void addConstituentPDF(PDF* pdf) {
      _pdfs.push_back(pdf);
    }

    /// Append several PDFs to the list of constituents
    template <typename PDFPTR>
    void addConstituentPDFs(std::vector<PDFPTR>& pdfs) {
      BOOST_FOREACH (PDF* p, pdfs) addConstituentPDF(p);
    }

    /// Clear the list of constituent PDFs, deleting the objects
    void reset() {
      BOOST_FOREACH (PDF* p, _pdfs) delete p;
      _pdfs.clear();
    }

    //@}


  protected:

    /// Get PDF xf(x,Q2) value (via multiplicative PDF combination)
    double _xfxQ2(int id, double x, double q2) const {
      double rtn = 1;
      BOOST_FOREACH (const PDF* p, constituentPDFs())
        rtn *= p->xfxQ2(id, x, q2);
      return rtn;
    }

    /// Collection of owned PDF pointers
    // std::vector< unique_ptr<PDF> > _pdfs;
    std::vector<PDF*> _pdfs;

  };


}
#endif
