// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

/// @file LHAGlue.h
/// A file that provides backwards compatibility for some C functions from LHAPDF 5.x


/// A special C++ function to return the PDF name + code currently being used via LHAGlue.
std::string lhaglue_get_current_pdf(int nset=1);



namespace LHAPDF {

  /// @brief Only Provided for LHAPDF5 compatibility. Distinction between evolution 
  /// or interpolation PDF sets.
  /// Enum to choose whether evolution (i.e. @c LHpdf data file) or
  /// interpolation (i.e. @c LHgrid data file) is used. 
  /// This distinction ismeaningless in LHAPDF6 interpolation is always used.
  enum SetType {
    EVOLVE = 0, LHPDF = 0,
    INTERPOLATE = 1, LHGRID = 1
  };

  /// check if the pdf uses photons
  bool hasPhoton();

  /// prepends path to path list
  void setPDFPath(const string& path);
  std::string pdfsetsPath();


  /// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDFSetByName(const std::string& filename);
  void initPDFSetByName(const std::string& filename, SetType type);

  /// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDFSetByName(int nset, const std::string& filename);
  void initPDFSetByName(int nset, const std::string& filename, SetType type);



  /// Number of members available in the current set.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int numberPDF();

  /// Number of members available in the current set.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int numberPDF(int nset);

  /// The choice of PDF member out of one distribution.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDF(int memset);

  /// The choice of PDF member out of one distribution.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDF(int nset, int memset);

  /// Nucleon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - flavour encoding as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double xfx(double x, double Q, int fl);

  /// Nucleon PDF: returns @c x f(x, Q) for flavour @a fl - flavour encoding as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double xfx(int nset, double x, double Q, int fl);

  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void xfx(double x, double Q, double* results);

  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void xfx(int nset, double x, double Q, double* results);

  /// Nucleon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  std::vector<double> xfx(double x, double Q);

  /// Nucleon PDF: returns a vector @c x f_i(x, Q) with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  std::vector<double> xfx(int nset, double x, double Q);




  /// Order of \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  int getOrderAlphaS();

  /// Number of flavours used in current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int getNf(int nset);


  /// Number of flavours used in current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  int getNf();


  /// Minimum X for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getXmin(int nmem);


  /// Minimum X for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getXmin(int nset, int nmem);

  /// Maximum X for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getXmax(int nset, int nmem);


  /// Maximum X for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getXmax(int nmem);

  /// Minimum Q2 for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getQ2min(int nset, int nmem);

  /// Minimum Q2 for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getQ2min(int nmem);

  /// Maximum Q2 for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getQ2max(int nset, int nmem);

  /// Maximum Q2 for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getQ2max(int nmem);


  /// Mass of quarks for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getQMass(int nset, int nf);

  /// Mass of quarks for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getQMass(int nf);


  /// Mass of quarks for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double getThreshold(int nset, int nf);

  /// Mass of quarks for current PDF
  /// @deprecated Use the proper C++ interface of LHAPDF6 -instead!
  double getThreshold(int nf);

  /// Order of \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  int getOrderAlphaS(int nset);


  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double alphasPDF(double Q);

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  double alphasPDF(int nset, double Q);

  /// @brief Use @a member in current PDF set.
  /// This operation is computationally cheap.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void usePDFMember(int member);

  /// @brief Use @a member in PDF set @a nset (multi-set version).
  /// This operation is computationally cheap.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void usePDFMember(int nset, int member);

  /// Initialise @a member in PDF set @a setid.
  /// @deprecated Use the proper C++ interface of LHAPDF6 instead!
  void initPDFSet(const std::string& name, int member=0);
  void initPDFSet(int nset ,const std::string& name, int nmem=0);
  void initPDFSet(int setid, int member=0);
  void initPDFSet(int nset , int setid, int nmem=0);
  /// Initialise @a member in PDF set @a name, of type @a type. - LHAPDF5 compatibility
  void initPDFSet(const std::string& name, SetType type, int member=0);  
  /// Initialise @a member in PDF set @a name, of type @a type (multi-set version) - LHAPDF5 compatibility
  void initPDFSet(int nset, const std::string& name, SetType type, int member=0);


}
