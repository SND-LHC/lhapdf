#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"

using namespace std;

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
    PDFSetHandler() : currentmem(1)
    { } //< It'll be stored in a map so we need one of these...

    /// Constructor from a PDF set name
    PDFSetHandler(const string& name)
      : setname(name)
    {
      loadMember(1);
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
      assert(mem > 0);
      if (members.find(mem) == members.end())
        members[mem] = PDFPtr(LHAPDF::mkPDF(setname, mem));
      currentmem = mem;
    }

    /// Actively delete a PDF member to save memory
    void unloadMember(int mem) {
      members.erase(mem);
      const int nextmem = (!members.empty()) ? members.begin()->first : 1;
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
      } catch (const std::exception& e) {
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
    } catch (const std::exception& e) {
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
    /// @todo Implement me! (and improve param names)
    throw LHAPDF::NotImplementedError("Photon structure functions are not yet supported");
  }

  /// Get xf(x) values for common partons from a photon PDF (non-multiset version)
  void evolvepdfp_(double& x, double& q, double& p2, int& ip2, double& fxq) {
    int nset1 = 1;
    evolvepdfpm_(nset1, x, q, p2, ip2, fxq);
  }


  /// @todo Doc
  void numberpdfm_(int& nset, int& numpdf) {
    /// @todo Implement me!
  }

  /// @todo Doc
  void numberpdf_(int& numpdf) {
    int nset1 = 1;
    numberpdfm_(nset1, numpdf);
  }


  /// @todo Doc
  void getnfm_(int& nset, double& nf) {
    /// @todo Implement me! (and improve param names)
  }

  /// @todo Doc
  void getnf_(double& nf){
    int nset1 = 1;
    getnfm_(nset1, nf);
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


  /// PDFLIB initialisation function
  void pdfset_(const char* par, const double* value, int parlength) {
    /// Take PDF ID from value[0]
    ACTIVESETS[1] = PDFSetHandler(value[0]);
    /// @todo How to use the par string?... most important for PYTHIA6?
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
