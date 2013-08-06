// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
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
      assert(mem >= 0);
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


  static const int LHAPIDS[13] = { -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 };
  static map<int, PDFSetHandler> ACTIVESETS;


}


string lhaglue_get_current_pdf(int nset) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    return "NONE";
  return ACTIVESETS[nset].activemember()->set().name() + " (" +
    LHAPDF::to_str(ACTIVESETS[nset].activemember()->lhapdfID()) + ")";
}


extern "C" {

  /// LHAPDF version
  void getlhapdfversion_(char *s, size_t len) {
    /// @todo Works? Need to check Fortran string return, string macro treatment, etc.
    strncpy(s, LHAPDF_VERSION, len);
  }


  /// Load a PDF set
  ///
  /// @todo Does this version actually take a *path*? What to do?
  void initpdfsetm_(int& nset, const char* setpath, int setpathlength) {
    // Strip file extension for backward compatibility
    /// @todo Don't we need to use substr & setnamelength since Fortran strs are not 0-terminated?
    const boost::filesystem::path p = setpath;
    string path = (p.extension().empty()) ? p.native() : p.stem().native(); //< @todo Will be wrong if a structured path is given
    /// Correct the misnamed CTEQ6L1/CTEQ6ll set name as a backward compatibility special case.
    if (boost::algorithm::to_lower_copy(path) == "cteq6ll") path = "CTEQ6L1";
    // Create the PDF set with index nset
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      ACTIVESETS[nset] = PDFSetHandler(path); //< @todo Will be wrong if a structured path is given
  }

  /// Load a PDF set (non-multiset version)
  void initpdfset_(const char* setpath, int setpathlength) {
    int nset1 = 1;
    initpdfsetm_(nset1, setpath, setpathlength);
  }


  /// Load a PDF set by name
  void initpdfsetbynamem_(int& nset, const char* setname, int setnamelength) {
    // Strip file extension for backward compatibility
    /// @todo Don't we need to use substr & setnamelength since Fortran strs are not 0-terminated?
    const boost::filesystem::path p = setname;
    string name = (p.extension().empty()) ? p.native() : p.stem().native();
    /// Correct the misnamed CTEQ6L1/CTEQ6ll set name as a backward compatibility special case.
    if (boost::algorithm::to_lower_copy(name) == "cteq6ll") name = "CTEQ6L1";
    // Create the PDF set with index nset
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      ACTIVESETS[nset] = PDFSetHandler(name);
  }

  /// Load a PDF set by name (non-multiset version)
  void initpdfsetbyname_(const char* setname, int setnamelength) {
    int nset1 = 1;
    initpdfsetbynamem_(nset1, setname, setnamelength);
  }


  /// Load a PDF in current set
  void initpdfm_(int& nset, int& nmember) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    ACTIVESETS[nset].loadMember(nmember);
  }

  /// Load a PDF in current set (non-multiset version)
  void initpdf_(int& nmember) {
    int nset1 = 1;
    initpdfm_(nset1, nmember);
  }


  /// Get xf(x) values for common partons from current PDF
  void evolvepdfm_(int& nset, double& x, double& q, double* fxq) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // Evaluate for the LHAPDF5 standard partons
    for (size_t i = 0; i < 13; ++i) {
      try {
        fxq[i] = ACTIVESETS[nset].activemember()->xfxQ(LHAPIDS[i], x, q);
      } catch (const exception& e) {
        fxq[i] = 0;
      }
    }
  }

  /// Get xf(x) values for common partons from current PDF (non-multiset version)
  void evolvepdf_(double& x, double& q, double* fxq) {
    int nset1 = 1;
    evolvepdfm_(nset1, x, q, fxq);
  }


  /// Determine if the current PDF has a photon flavour (historically only MRST2004QED)
  /// @todo Function rather than subroutine?
  bool has_photon_() {
    /// @todo Only apply to nset = 1? Or do we need to track the current nset value?
    return ACTIVESETS[1].activemember()->hasFlavor(22);
  }


  /// Get xfx values from current PDF, including an extra photon flavour
  void evolvepdfphotonm_(int& nset, double& x, double& q, double* fxq, double& photonfxq) {
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
  }

  /// Get xfx values from current PDF, including an extra photon flavour (non-multiset version)
  void evolvepdfphoton_(double& x, double& q, double* fxq, double& photonfxq) {
    int nset1 = 1;
    evolvepdfphotonm_(nset1, x, q, fxq, photonfxq);
  }


  /// Get xf(x) values for common partons from a photon PDF
  void evolvepdfpm_(int& nset, double& x, double& q, double& p2, int& ip2, double& fxq) {
    /// @todo Implement me!
    throw LHAPDF::NotImplementedError("Photon structure functions are not yet supported");
  }

  /// Get xf(x) values for common partons from a photon PDF (non-multiset version)
  void evolvepdfp_(double& x, double& q, double& p2, int& ip2, double& fxq) {
    int nset1 = 1;
    evolvepdfpm_(nset1, x, q, p2, ip2, fxq);
  }


  void numberpdfm_(int& nset, int& numpdf) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // set equal to the number of members  for the requested set
    numpdf=  ACTIVESETS[nset].activemember()->info().get_entry_as<int>("NumMembers");
    // reproduce old lhapdf v5 behaviour , subtract 1 if more than 1 member set
    if (numpdf > 1) numpdf-=1;
  }

  void numberpdf_(int& numpdf) {
    int nset1 = 1;
    numberpdfm_(nset1, numpdf);
  }


  /// Get the number of flavours
  void getnfm_(int& nset, double& nf) {
    /// @todo Implement me! (and improve param names)
  }

  /// Get the number of flavours (non-multiset version)
  void getnf_(double& nf) {
    int nset1 = 1;
    getnfm_(nset1, nf);
  }


  /// @todo Doc and better arg name
  void lhaprint_(int& a){
    // Do nothing for now
    /// @todo Can this be mapped?
  }


  /// @todo Need getnset_ and getnmem_


  /// @brief Set LHAPDF parameters
  ///
  /// Does nothing for now
  void setlhaparm_(const char* par, int parlength) {
    /// @todo Can any Fortran LHAPDF params be usefully mapped?
  }


  /// PDFLIB initialisation function
  void pdfset_(const char* par, const double* value, int parlength) {

    // Initialise struct equivalents to common blocks with sensible values.
    w50512_.qcdl4 = 0.215;
    w50512_.qcdl5 = 0.165;
    w50513_.xmin = 0.0;
    w50513_.xmax = 1.0;
    w50513_.q2min = 1.0;
    w50513_.q2max = 1.0e5;
    lhapdfr_.qcdlha4 = 0.0;
    lhapdfr_.qcdlha5 = 0.0;
    lhapdfr_.nfllha = 4;

    string my_par(par);
    if (my_par.find("NPTYPE") != string::npos) {
      /// @todo Remove noisiness? Useful for now
      cout << "==== PYTHIA WILL USE LHAPDFv6 ====" << endl;
    } else if (my_par.find("HWLHAPDF") != string::npos) {
      /// @todo Remove noisiness? Useful for now
      cout << "==== HERWIG WILL USE LHAPDFv6 ====" << endl;
    }

    // Take PDF ID from value[2]
    ACTIVESETS[1] = PDFSetHandler(value[2]+1000*value[1]);

    // Need to extract parameters for common blocks
    PDFPtr pdf = ACTIVESETS[1].activemember();

    w50513_.xmin=pdf->info().get_entry_as<double>("XMin");
    w50513_.xmax=pdf->info().get_entry_as<double>("XMax");
    w50513_.q2min=pdf->info().get_entry_as<double>("Q2Min");
    w50513_.q2max=pdf->info().get_entry_as<double>("Q2Max");
    w50512_.qcdl4=pdf->info().get_entry_as<double>("AlphaS_Lambda4");
    w50512_.qcdl5=pdf->info().get_entry_as<double>("AlphaS_Lambda5");
    lhapdfr_.qcdlha4=pdf->info().get_entry_as<double>("AlphaS_Lambda4");
    lhapdfr_.qcdlha5=pdf->info().get_entry_as<double>("AlphaS_Lambda5");
    // BEGIN: used to test behaviour versus lhapdf 5.x
    //    w50512_.qcdl4=0.192;
    //    w50512_.qcdl5=0.192;
    //    lhapdfr_.qcdlha4=0.192;
    //    lhapdfr_.qcdlha5=0.192;
    // END:  backwards compatibility test
  }


  /// PDFLIB nucleon structure function querying
  void structm_(const double& x, const double& q,
                double& upv, double& dnv, double& usea, double& dsea,
                double& str, double& chm, double& bot, double& top, double& glu) {
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
    /// @todo Do anything?
  }


}


// LHAPDF namespace compatibility code

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


void LHAPDF::initPDFSet(const string& filename) {
  initPDFSetByName(filename);
}



#define SIZE 999
void initPDFSetByName(const string& filename) {
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbyname_(cfilename, filename.length());
}

void initPDFSetByName(int nset, const string& filename) {
  char cfilename[SIZE+1];
  strncpy(cfilename, filename.c_str(), SIZE);
  initpdfsetbynamem_(nset, cfilename, filename.length());
}


double LHAPDF::alphasPDF(double Q) {
  return LHAPDF::alphasPDF(1, Q) ;
}

double LHAPDF::alphasPDF(int nset, double Q) {
  if (ACTIVESETS.find(nset) == ACTIVESETS.end())
    throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
  // return alphaS for the requested set
  return ACTIVESETS[nset].activemember()->alphasQ(Q);
}

void LHAPDF::usePDFMember(int member) {
  initpdf_(member);
}

void LHAPDF::usePDFMember(int nset, int member) {
  initpdfm_(nset, member);
}
