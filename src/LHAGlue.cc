// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"
#include "LHAPDF/Version.h"
#include "LHAPDF/LHAGlue.h"

using namespace std;


// We have to create and initialise some common blocks here for backwards compatibility
struct w50512 {
  double qcdl4, qcdl5;
};
w50512 w50512_;

struct w50513 {
  double xmin, xmax, q2min, q2max;
};
w50513 w50513_;

struct lhapdfr {
  double qcdlha4, qcdlha5;
  int nfllha;
};
lhapdfr lhapdfr_;



namespace { //< Unnamed namespace to restrict visibility to this file

  /// @brief PDF object storage here is a smart pointer to ensure deletion of created PDFs
  ///
  /// NB. std::auto_ptr cannot be stored in STL containers, hence we use
  /// boost::shared_ptr. std::unique_ptr is the nature replacement when C++11
  /// is globally available.
  typedef boost::shared_ptr<LHAPDF::PDF> PDFPtr;

  /// @brief A struct for handling the active PDFs for the Fortran interface.
  ///
  /// We operate in a string-based way, since maybe there will be sets with names, but no
  /// index entry in pdfsets.index.
  ///
  /// @todo Can we avoid the strings and just work via the LHAPDF ID and factory construction?
  ///
  /// Smart pointers are used in the native map used for PDF member storage so
  /// that they auto-delete if the PDFSetHandler that holds them goes out of
  /// scope (i.e. is overwritten).
  struct PDFSetHandler {

    /// Default constructor
    PDFSetHandler() : currentmem(0)
    { } //< It'll be stored in a map so we need one of these...

    /// Constructor from a PDF set name
    PDFSetHandler(const string& name)
      : setname(name)
    {
      loadMember(0);
    }

    /// Constructor from a PDF set's LHAPDF ID code
    PDFSetHandler(int lhaid) {
      pair<string,int> set_mem = LHAPDF::lookupPDF(lhaid);
      setname = set_mem.first;
      loadMember(set_mem.second);
    }

    /// @brief Load a new PDF member
    ///
    /// If it's already loaded, the existing object will not be reloaded.
    void loadMember(int mem) {
      if (mem < 0)
        throw LHAPDF::UserError("Tried to load a negative PDF member ID: " + LHAPDF::to_str(mem) + " in set " + setname);
      if (members.find(mem) == members.end())
        members[mem] = PDFPtr(LHAPDF::mkPDF(setname, mem));
      currentmem = mem;
    }

    /// Actively delete a PDF member to save memory
    void unloadMember(int mem) {
      members.erase(mem);
      const int nextmem = (!members.empty()) ? members.begin()->first : 0;
      loadMember(nextmem);
    }

    /// @brief Get a PDF member
    ///
    /// Non-const because it can secretly load the member. Not that constness
    /// matters in a Fortran interface utility function!
    const PDFPtr member(int mem) {
      loadMember(mem);
      return members.find(mem)->second;
    }

    /// Get the currently active PDF member
    ///
    /// Non-const because it can secretly load the member. Not that constness
    /// matters in a Fortran interface utility function!
    const PDFPtr activemember() {
      return member(currentmem);
    }

    /// The currently active member in this set
    int currentmem;

    /// Name of this set
    string setname;

    /// Map of pointers to selected member PDFs
    ///
    // /// It's mutable so that a "const" member-getting operation can implicitly
    // /// load a new PDF object. Good idea / bad idea? Disabled for now.
    // mutable map<int, PDFPtr> members;
    map<int, PDFPtr> members;
  };


  /// Collection of active sets
  static map<int, PDFSetHandler> ACTIVESETS;

  /// The currently active set
  int currentset = 0;

}



string lhaglue_get_current_pdf(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    return "NONE";
  currentset = nset;
  return ACTIVESETS[nset].activemember()->set().name() + " (" +
    LHAPDF::to_str(ACTIVESETS[nset].activemember()->lhapdfID()) + ")";
}



extern "C" {


  /// LHAPDF version
  void getlhapdfversion_(char* s, size_t len) {
    /// @todo Works? Need to check Fortran string return, string macro treatment, etc.
    strncpy(s, LHAPDF_VERSION, len);
  }


  /// Return a dummy max number of sets (there is no limitation now)
  void getmaxnumsets_(int& nmax) {
    nmax = 1000;
  }


  /// Set path
  void setpdfpath_(const char* s, size_t len) {
    /// @todo Works? Need to check Fortran string return, C-string copying, null termination
    char s2[1024];
    s2[len] = '\0';
    strncpy(s2, s, len);
    LHAPDF::pathsPrepend(s2);
  }


  /// Load a PDF set
  ///
  /// @todo Does this version actually take a *path*? What to do?
  void initpdfsetm_(const int& nset, const char* setpath, int setpathlength) {
    // Strip file extension for backward compatibility
    string fullp = string(setpath, setpathlength);
    // Remove trailing whitespace
    fullp.erase( std::remove_if( fullp.begin(), fullp.end(), ::isspace ), fullp.end() );
    // Use only items after the last /
    const string pap = LHAPDF::dirname(fullp);
    const string p = LHAPDF::basename(fullp);
    // Prepend path to search area
    LHAPDF::pathsPrepend(pap);
    // Handle extensions
    string path = LHAPDF::file_extn(p).empty() ? p : LHAPDF::file_stem(p);
    /// @note We correct the misnamed CTEQ6L1/CTEQ6ll set name as a backward compatibility special case.
    if (boost::algorithm::to_lower_copy(path) == "cteq6ll") path = "cteq6l1";
    // Create the PDF set with index nset
    // if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    ACTIVESETS[nset] = PDFSetHandler(path); //< @todo Will be wrong if a structured path is given
    currentset = nset;
  }
  /// Load a PDF set (non-multiset version)
  void initpdfset_(const char* setpath, int setpathlength) {
    int nset1 = 1;
    initpdfsetm_(nset1, setpath, setpathlength);
  }


  /// Load a PDF set by name
  void initpdfsetbynamem_(const int& nset, const char* setname, int setnamelength) {
    // Truncate input to size setnamelength
    string p = setname;
    p.erase(setnamelength, std::string::npos);
    // Strip file extension for backward compatibility
    string name = LHAPDF::file_extn(p).empty() ? p : LHAPDF::file_stem(p);
    // Remove trailing whitespace
    name.erase( std::remove_if( name.begin(), name.end(), ::isspace ), name.end() );
    /// @note We correct the misnamed CTEQ6L1/CTEQ6ll set name as a backward compatibility special case.
    if (boost::algorithm::to_lower_copy(name) == "cteq6ll") name = "cteq6l1";
    // Create the PDF set with index nset
    // if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    ACTIVESETS[nset] = PDFSetHandler(name);
    currentset = nset;
  }
  /// Load a PDF set by name (non-multiset version)
  void initpdfsetbyname_(const char* setname, int setnamelength) {
    int nset1 = 1;
    initpdfsetbynamem_(nset1, setname, setnamelength);
  }


  /// Load a PDF in current set
  void initpdfm_(const int& nset, const int& nmember) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    ACTIVESETS[nset].loadMember(nmember);
    currentset = nset;
  }
  /// Load a PDF in current set (non-multiset version)
  void initpdf_(const int& nmember) {
    int nset1 = 1;
    initpdfm_(nset1, nmember);
  }


  void getnset_(int& nset) {
    nset = currentset;
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  }

  void setnset_(const int& nset) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" +
                              LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
  }

  void getnmem_(int& nset, int& nmem) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    nmem = ACTIVESETS[nset].currentmem;
  }

  void setnmem_(const int& nset, const int& nmem) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" +
                              LHAPDF::to_str(nset) + " but it is not initialised");
    ACTIVESETS[nset].loadMember(nmem);
  }


  /// Get xf(x) values for common partons from current PDF
  void evolvepdfm_(const int& nset, const double& x, const double& q, double* fxq) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // Evaluate for the 13 LHAPDF5 standard partons (-6..6)
    for (size_t i = 0; i < 13; ++i) {
      try {
        fxq[i] = ACTIVESETS[nset].activemember()->xfxQ(i-6, x, q);
      } catch (const exception& e) {
        fxq[i] = 0;
      }
    }
  }
  /// Get xf(x) values for common partons from current PDF (non-multiset version)
  void evolvepdf_(const double& x, const double& q, double* fxq) {
    int nset1 = 1;
    evolvepdfm_(nset1, x, q, fxq);
  }


  /// Determine if the current PDF has a photon flavour (historically only MRST2004QED)
  /// @todo Function rather than subroutine?
  bool has_photon_() {
    /// @todo Only apply to nset = 1? Or do we need to track the current nset value?
    currentset = 1;
    return ACTIVESETS[1].activemember()->hasFlavor(22);
  }


  /// Get xfx values from current PDF, including an extra photon flavour
  void evolvepdfphotonm_(const int& nset, const double& x, const double& q, double* fxq, double& photonfxq) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // First evaluate the "normal" partons
    evolvepdfm_(nset, x, q, fxq);
    // Then evaluate the photon flavor (historically only for MRST2004QED)
    try {
      photonfxq = ACTIVESETS[nset].activemember()->xfxQ(22, x, q);
    } catch (const exception& e) {
      photonfxq = 0;
    }
    currentset = nset;
  }
  /// Get xfx values from current PDF, including an extra photon flavour (non-multiset version)
  void evolvepdfphoton_(const double& x, const double& q, double* fxq, double& photonfxq) {
    int nset1 = 1;
    evolvepdfphotonm_(nset1, x, q, fxq, photonfxq);
  }


  /// Get xf(x) values for common partons from a photon PDF
  void evolvepdfpm_(const int& nset, const double& x, const double& q, const double& p2, const int& ip2, double& fxq) {
    /// @todo Implement me!
    throw LHAPDF::NotImplementedError("Photon structure functions are not yet supported");
  }
  /// Get xf(x) values for common partons from a photon PDF (non-multiset version)
  void evolvepdfp_(const double& x, const double& q, const double& p2, const int& ip2, double& fxq) {
    int nset1 = 1;
    evolvepdfpm_(nset1, x, q, p2, ip2, fxq);
  }


  /// Get the number of members in the set
  void numberpdfm_(const int& nset, int& numpdf) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    // set equal to the number of members  for the requested set
    numpdf=  ACTIVESETS[nset].activemember()->info().get_entry_as<int>("NumMembers");
    // reproduce old lhapdf v5 behaviour , subtract 1 if more than 1 member set
    if (numpdf > 1) numpdf-=1;
  }
  /// Get the number of members in the set (non-multiset version)
  void numberpdf_(int& numpdf) {
    int nset1 = 1;
    numberpdfm_(nset1, numpdf);
  }


  /// Get the alpha_s order for the set
  void getorderasm_(const int& nset, int& oas) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    // Set equal to the number of members for the requested set
    oas = ACTIVESETS[nset].activemember()->info().get_entry_as<int>("AlphaS_OrderQCD");
  }
  /// Get the alpha_s order for the set (non-multiset version)
  void getorderas_(int& oas) {
    int nset1 = 1;
    getorderasm_(nset1, oas);
  }


  /// Get the number of flavours
  void getnfm_(const int& nset, int& nf) {
    //nf = ACTIVESETS[nset].activemember()->info().get_entry_as<int>("AlphaS_NumFlavors");
    nf = ACTIVESETS[nset].activemember()->info().get_entry_as<int>("NumFlavors");
    currentset = nset;
  }
  /// Get the number of flavours (non-multiset version)
  void getnf_(int& nf) {
    int nset1 = 1;
    getnfm_(nset1, nf);
  }


  /// Get nf'th quark mass
  void getqmassm_(const int& nset, const int& nf, double& mass) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    if (nf*nf == 1) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MDown");
    else if (nf*nf == 4) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MUp");
    else if (nf*nf == 9) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MStrange");
    else if (nf*nf == 16) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MCharm");
    else if (nf*nf == 25) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MBottom");
    else if (nf*nf == 36) mass = ACTIVESETS[nset].activemember()->info().get_entry_as<double>("MTop");
    else throw LHAPDF::UserError("Trying to get quark mass for invalid quark ID #" + LHAPDF::to_str(nf));
  }
  /// Get nf'th quark mass (non-multiset version)
  void getqmass_(const int& nf, double& mass) {
    int nset1 = 1;
    getqmassm_(nset1, nf, mass);
  }


  /// Get the nf'th quark threshold
  /// @todo Update when we can support generalized thresholds
  void getthresholdm_(const int& nset, const int& nf, double& Q) {
    currentset = nset;
    getqmassm_(nset, nf, Q);
  }
  /// Get the nf'th quark threshold
  void getthreshold_(const int& nf, double& Q) {
    int nset1 = 1;
    getthresholdm_(nset1, nf, Q);
  }


  /// @todo Doc and better arg name
  void lhaprint_(int& a){
    // Do nothing for now
    /// @todo Can this be mapped?
  }


  /// @brief Set LHAPDF parameters
  ///
  /// Does nothing for now
  void setlhaparm_(const char* par, int parlength) {
    /// @todo Can any Fortran LHAPDF params be usefully mapped?
  }


  double alphaspdfm_(const int& nset, const double& Q){
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    return ACTIVESETS[nset].activemember()->alphasQ(Q);
  }
  double alphaspdf_(const double& Q){
    int nset1 = 1;
    return alphaspdfm_(nset1, Q);
  }


  /// PDFLIB initialisation function
  void pdfset_(const char* par, const double* value, int parlength) {

    // Identify the calling program (yuck!)
    string my_par(par);
    if (my_par.find("NPTYPE") != string::npos) {
      /// @todo Remove noisiness? Useful for now
      cout << "==== PYTHIA WILL USE LHAPDFv6 ====" << endl;
      // Take PDF ID from value[2]
      ACTIVESETS[1] = PDFSetHandler(value[2]+1000*value[1]);
    } else if (my_par.find("HWLHAPDF") != string::npos) {
      /// @todo Remove noisiness? Useful for now
      cout << "==== HERWIG WILL USE LHAPDFv6 ====" << endl;
      // Take PDF ID from value[0]
      ACTIVESETS[1] = PDFSetHandler(value[0]);
    } else if (my_par.find("DEFAULT") != string::npos) {
      cout << "==== POWHEG (OR POSSIBLY HERWIG) WILL USE LHAPDFv6 ====" << endl;
      // Take PDF ID from value[0]
      ACTIVESETS[1] = PDFSetHandler(value[0]);
    } else {
      cout << "==== UNKNOWN GENERATOR WILL USE LHAPDFv6 IN PDFLIB STYLE ====" << endl;
      // Take PDF ID from value[2]
      ACTIVESETS[1] = PDFSetHandler(value[2]+1000*value[1]);
    }

    currentset = 1;

    // Extract parameters for common blocks (with sensible fallback values)
    PDFPtr pdf = ACTIVESETS[1].activemember();
    w50513_.xmin = pdf->info().get_entry_as<double>("XMin", 0.0);
    w50513_.xmax = pdf->info().get_entry_as<double>("XMax", 1.0);
    w50513_.q2min = LHAPDF::sqr(pdf->info().get_entry_as<double>("QMin", 1.0));
    w50513_.q2max = LHAPDF::sqr(pdf->info().get_entry_as<double>("QMax", 1.0e5));
    w50512_.qcdl4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4", 0.0);
    w50512_.qcdl5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5", 0.0);
    lhapdfr_.qcdlha4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4", 0.0);
    lhapdfr_.qcdlha5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5", 0.0);
    lhapdfr_.nfllha = 4;
    // Activate legacy/compatibility LHAPDF5-type behaviour re. broken Lambda values
    if (pdf->info().get_entry_as<bool>("Pythia6LambdaV5Compat", true)) {
      w50512_.qcdl4 = 0.192;
      w50512_.qcdl5 = 0.192;
      lhapdfr_.qcdlha4 = 0.192;
      lhapdfr_.qcdlha5 = 0.192;
    }
  }

  /// PDFLIB nucleon structure function querying
  void structm_(const double& x, const double& q,
                double& upv, double& dnv, double& usea, double& dsea,
                double& str, double& chm, double& bot, double& top, double& glu) {
    currentset = 1;
    /// Fill (partial) parton return variables
    PDFPtr pdf = ACTIVESETS[1].activemember();
    dsea = pdf->xfxQ(-1, x, q);
    usea = pdf->xfxQ(-2, x, q);
    dnv = pdf->xfxQ(1, x, q) - dsea;
    upv = pdf->xfxQ(2, x, q) - usea;
    str = pdf->xfxQ(3, x, q);
    chm = (pdf->hasFlavor(4)) ? pdf->xfxQ(4, x, q) : 0;
    bot = (pdf->hasFlavor(5)) ? pdf->xfxQ(5, x, q) : 0;
    top = (pdf->hasFlavor(6)) ? pdf->xfxQ(6, x, q) : 0;
    glu = pdf->xfxQ(21, x, q);
  }

  /// PDFLIB photon structure function querying
  void structp_(const double& x, const double& q2, const double& p2, const double& ip2,
                double& upv, double& dnv, double& usea, double& dsea,
                double& str, double& chm, double& bot, double& top, double& glu) {
    throw LHAPDF::NotImplementedError("Photon structure functions are not yet supported");
  }

  /// PDFLIB statistics on PDF under/overflows
  void pdfsta_() {
    /// @note Can't do anything...
  }


  void getxmin_(const int& nmem, double& xmin) {
    currentset = 1;
    ACTIVESETS[1].loadMember(nmem);
    xmin = ACTIVESETS[1].activemember()->info().get_entry_as<double>("XMin");
  }

  void getxmax_(const int& nmem, double& xmax) {
    currentset = 1;
    ACTIVESETS[1].loadMember(nmem);
    xmax = ACTIVESETS[1].activemember()->info().get_entry_as<double>("XMax");
  }

  void getq2min_(const int& nmem, double& q2min) {
    currentset = 1;
    ACTIVESETS[1].loadMember(nmem);
    q2min = LHAPDF::sqr(ACTIVESETS[1].activemember()->info().get_entry_as<double>("QMin"));
  }

  void getq2max_(const int& nmem, double& q2max) {
    currentset = 1;
    ACTIVESETS[1].loadMember(nmem);
    q2max = LHAPDF::sqr(ACTIVESETS[1].activemember()->info().get_entry_as<double>("QMax"));
  }


  /// Backwards compatibility functions for LHAPDF5 calculations of
  /// PDF uncertainties and PDF correlations (G. Watt, March 2014).

  // subroutine GetPDFUncTypeM(nset,lMonteCarlo,lSymmetric)
  void getpdfunctypem_(const int& nset, int& lmontecarlo, int& lsymmetric) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    const string errorType = ACTIVESETS[nset].activemember()->set().errorType();
    if (errorType == "replicas") { // Monte Carlo PDF sets
      lmontecarlo = 1;
      lsymmetric = 1;
    } else if (errorType == "symmhessian") { // symmetric eigenvector PDF sets
      lmontecarlo = 0;
      lsymmetric = 1;
    } else { // default: assume asymmetric Hessian eigenvector PDF sets
      lmontecarlo = 0;
      lsymmetric = 0;
    }
  }

  // subroutine GetPDFUncType(lMonteCarlo,lSymmetric)
  void getpdfunctype_(int& lmontecarlo, int& lsymmetric) {
    int nset1 = 1;
    getpdfunctypem_(nset1, lmontecarlo, lsymmetric);
  }

  // subroutine GetPDFuncertaintyM(nset,values,central,errplus,errminus,errsym)
  void getpdfuncertaintym_(const int& nset, const double* values, double& central, double& errplus, double& errminus, double& errsymm) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    const size_t nmem = ACTIVESETS[nset].activemember()->set().size()-1;
    const vector<double> vecvalues(values, values + nmem + 1);
    LHAPDF::PDFUncertainty err = ACTIVESETS[nset].activemember()->set().uncertainty(vecvalues, -1);
    central = err.central;
    errplus = err.errplus;
    errminus = err.errminus;
    errsymm = err.errsymm;
  }

  // subroutine GetPDFuncertainty(values,central,errplus,errminus,errsym)
  void getpdfuncertainty_(const double* values, double& central, double& errplus, double& errminus, double& errsymm) {
    int nset1 = 1;
    getpdfuncertaintym_(nset1, values, central, errplus, errminus, errsymm);
  }

  // subroutine GetPDFcorrelationM(nset,valuesA,valuesB,correlation)
  void getpdfcorrelationm_(const int& nset, const double* valuesA, const double* valuesB, double& correlation) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    currentset = nset;
    const size_t nmem = ACTIVESETS[nset].activemember()->set().size()-1;
    const vector<double> vecvaluesA(valuesA, valuesA + nmem + 1);
    const vector<double> vecvaluesB(valuesB, valuesB + nmem + 1);
    correlation = ACTIVESETS[nset].activemember()->set().correlation(vecvaluesA,vecvaluesB);
  }

  // subroutine GetPDFcorrelation(valuesA,valuesB,correlation)
  void getpdfcorrelation_(const double* valuesA, const double* valuesB, double& correlation) {
    int nset1 = 1;
    getpdfcorrelationm_(nset1, valuesA, valuesB, correlation);
  }


}


// LHAPDF namespace C++ compatibility code
#ifdef ENABLE_LHAGLUE_CXX


void LHAPDF::setVerbosity(LHAPDF::Verbosity noiselevel) {
  LHAPDF::setVerbosity((int) noiselevel);
}

void LHAPDF::setPDFPath(const string& path) {
  pathsPrepend(path);
}

string LHAPDF::pdfsetsPath() {
  return paths()[0];
}

int LHAPDF::numberPDF() {
  int nmem;
  numberpdf_(nmem);
  return nmem;
}
int LHAPDF::numberPDF(int nset) {
  int nmem;
  numberpdfm_(nset,nmem);
  return nmem;
}

void LHAPDF::initPDF( int memset) {
  int nset1 = 1;
  initpdfm_(nset1, memset);
}
void LHAPDF::initPDF(int nset, int memset) {
  initpdfm_(nset, memset);
}


double LHAPDF::xfx(double x, double Q, int fl) {
  vector<double> r(13);
  evolvepdf_(x, Q, &r[0]);
  return r[fl+6];
}
double LHAPDF::xfx(int nset, double x, double Q, int fl) {
  vector<double> r(13);
  evolvepdfm_(nset, x, Q, &r[0]);
  return r[fl+6];
}

vector<double> LHAPDF::xfx(double x, double Q) {
  vector<double> r(13);
  evolvepdf_(x, Q, &r[0]);
  return r;
}
vector<double> LHAPDF::xfx(int nset, double x, double Q) {
  vector<double> r(13);
  evolvepdfm_(nset, x, Q, &r[0]);
  return r;
}

void LHAPDF::xfx(double x, double Q, double* results) {
  evolvepdf_(x, Q, results);
}
void LHAPDF::xfx(int nset, double x, double Q, double* results) {
  evolvepdfm_(nset, x, Q, results);
}


vector<double> LHAPDF::xfxphoton(double x, double Q) {
  vector<double> r(13);
  double mphoton;
  evolvepdfphoton_(x, Q, &r[0], mphoton);
  r.push_back(mphoton);
  return r;
}
vector<double> LHAPDF::xfxphoton(int nset, double x, double Q) {
  vector<double> r(13);
  double mphoton;
  evolvepdfphotonm_(nset, x, Q, &r[0], mphoton);
  r.push_back(mphoton);
  return r;
}

void LHAPDF::xfxphoton(double x, double Q, double* results) {
  evolvepdfphoton_(x, Q, results, results[13]);
}
void LHAPDF::xfxphoton(int nset, double x, double Q, double* results) {
  evolvepdfphotonm_(nset, x, Q, results, results[13]);
}

double LHAPDF::xfxphoton(double x, double Q, int fl) {
  vector<double> r(13);
  double mphoton;
  evolvepdfphoton_(x, Q, &r[0], mphoton);
  if (fl == 7) return mphoton;
  return r[fl+6];
}
double LHAPDF::xfxphoton(int nset, double x, double Q, int fl) {
  vector<double> r(13);
  double mphoton;
  evolvepdfphotonm_(nset, x, Q, &r[0], mphoton);
  if ( fl == 7 ) return mphoton;
  return r[fl+6];
}


void LHAPDF::initPDFSet(const string& filename, int nmem) {
  initPDFSet(1,filename, nmem);
}

void LHAPDF::initPDFSet(int nset, const string& filename, int nmem) {
  initPDFSetByName(nset,filename);
  ACTIVESETS[nset].loadMember(nmem);
  currentset = nset;
}


void LHAPDF::initPDFSet(const string& filename, SetType type ,int nmem) {
  // silently ignore type
  initPDFSet(1,filename, nmem);
}

void LHAPDF::initPDFSet(int nset, const string& filename, SetType type ,int nmem) {
  // silently ignore type
  initPDFSetByName(nset,filename);
  ACTIVESETS[nset].loadMember(nmem);
  currentset = nset;
}

void LHAPDF::initPDFSet(int nset, int setid, int nmem) {
  ACTIVESETS[nset] = PDFSetHandler(setid); //
  currentset = nset;
}

void LHAPDF::initPDFSet(int setid, int nmem) {
  initPDFSet(1,setid,nmem);
}

#define SIZE 999
void LHAPDF::initPDFSetByName(const string& filename) {
  std::cout << "initPDFSetByName: " << filename << std::endl;
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbyname_(cfilename, filename.length());
}

void LHAPDF::initPDFSetByName(int nset, const string& filename) {
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbynamem_(nset, cfilename, filename.length());
}

void LHAPDF::initPDFSetByName(const string& filename, SetType type) {
  //silently ignore type
  std::cout << "initPDFSetByName: " << filename << std::endl;
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbyname_(cfilename, filename.length());
}

void LHAPDF::initPDFSetByName(int nset, const string& filename, SetType type) {
  //silently ignore type
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbynamem_(nset, cfilename, filename.length());
}


void LHAPDF::getDescription() {
  getDescription(1);
}

void LHAPDF::getDescription(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  cout << ACTIVESETS[nset].activemember()->set().description() << endl;
}


double LHAPDF::alphasPDF(double Q) {
  return LHAPDF::alphasPDF(1, Q) ;
}

double LHAPDF::alphasPDF(int nset, double Q) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS for the requested set
  return ACTIVESETS[nset].activemember()->alphasQ(Q);
}


bool LHAPDF::hasPhoton(){
  return has_photon_();
}


int LHAPDF::getOrderAlphaS() {
  return LHAPDF::getOrderAlphaS(1) ;
}

int LHAPDF::getOrderAlphaS(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  return ACTIVESETS[nset].activemember()->info().get_entry_as<int>("AlphaS_OrderQCD", -1);
}


int LHAPDF::getOrderPDF() {
  return LHAPDF::getOrderPDF(1) ;
}

int LHAPDF::getOrderPDF(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return PDF order for the requested set
  return ACTIVESETS[nset].activemember()->info().get_entry_as<int>("OrderQCD", -1);
}


double LHAPDF::getLam4(int nmem) {
  return LHAPDF::getLam4(1, nmem) ;
}

double LHAPDF::getLam4(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  ACTIVESETS[nset].loadMember(nmem);
  return ACTIVESETS[nset].activemember()->info().get_entry_as<double>("AlphaS_Lambda4", -1.0);
}


double LHAPDF::getLam5(int nmem) {
  return LHAPDF::getLam5(1, nmem) ;
}

double LHAPDF::getLam5(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  ACTIVESETS[nset].loadMember(nmem);
  return ACTIVESETS[nset].activemember()->info().get_entry_as<double>("AlphaS_Lambda5", -1.0);
}


int LHAPDF::getNf() {
  return LHAPDF::getNf(1) ;
}

int LHAPDF::getNf(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  return ACTIVESETS[nset].activemember()->info().get_entry_as<int>("NumFlavors");
}


double LHAPDF::getXmin(int nmem) {
  return LHAPDF::getXmin(1, nmem) ;
}

double LHAPDF::getXmin(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  ACTIVESETS[nset].loadMember(nmem);
  return ACTIVESETS[nset].activemember()->info().get_entry_as<double>("XMin");
}

double LHAPDF::getXmax(int nmem) {
  return LHAPDF::getXmax(1, nmem) ;
}

double LHAPDF::getXmax(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  ACTIVESETS[nset].loadMember(nmem);
  return ACTIVESETS[nset].activemember()->info().get_entry_as<double>("XMax");
}

double LHAPDF::getQ2min(int nmem) {
  return LHAPDF::getQ2min(1, nmem) ;
}

double LHAPDF::getQ2min(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  ACTIVESETS[nset].loadMember(nmem);
  return pow(ACTIVESETS[nset].activemember()->info().get_entry_as<double>("QMin"),2);
}

double LHAPDF::getQ2max(int nmem) {
  return LHAPDF::getQ2max(1,nmem) ;
}

double LHAPDF::getQ2max(int nset, int nmem) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  currentset = nset;
  // return alphaS Order for the requested set
  ACTIVESETS[nset].loadMember(nmem);
  return pow(ACTIVESETS[nset].activemember()->info().get_entry_as<double>("QMax"),2);
}

double LHAPDF::getQMass(int nf) {
  return LHAPDF::getQMass(1, nf) ;
}

double LHAPDF::getQMass(int nset, int nf) {
  double mass;
  getqmassm_(nset, nf, mass);
  return mass;
}

double LHAPDF::getThreshold(int nf) {
  return LHAPDF::getThreshold(1, nf) ;
}

double LHAPDF::getThreshold(int nset, int nf) {
  currentset = nset;
  /// @todo Update when we can support generalized thresholds
  return LHAPDF::getQMass(nset, nf) ;
}

void LHAPDF::usePDFMember(int member) {
  initpdf_(member);
}

void LHAPDF::usePDFMember(int nset, int member) {
  initpdfm_(nset, member);
}

#endif // ENABLE_LHAGLUE_CXX
