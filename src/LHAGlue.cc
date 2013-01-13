#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/Utils.h"

using namespace std;

namespace { //< Unnamed namespace to restrict visibility to this file

  /// A struct for handling the active PDFs for the Fortran interface.
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

    /// Internal storage is a smart pointer to ensure deletion of created PDFs
    ///
    /// NB. std::auto_ptr cannot be stored in STL containers, hence we use
    /// boost::shared_ptr. std::unique_ptr is the nature replacement when C++11
    /// is globally available.
    typedef boost::shared_ptr<LHAPDF::PDF> PDFPtr;

    /// Default constructor
    PDFSetHandler() { } //< It is stored in a map so we need one of these...

    /// Constructor from a PDF set name
    PDFSetHandler(const string& name)
      : currentmem(-1), setname(setname)
    {    }

    /// Constructor from a PDF set's LHAPDF ID code
    PDFSetHandler(int lhaid) {
      pair<string,int> set_mem = LHAPDF::lookupPDF(lhaid);
      setname = set_mem.first;
      loadMember(set_mem.second);
    }

    /// Load a new PDF member
    ///
    /// If it's already loaded, the existing object will be deleted and
    /// recreated. Due to the underlying auto_ptr mechanism no memory should be
    /// leaked as a result.
    void loadMember(int mem) {
      members[mem] = PDFPtr(LHAPDF::mkPDF(setname, mem));
      currentmem = mem;
    }

    /// Actively delete a PDF member to save memory
    void unloadMember(int mem) {
      members.erase(mem);
      currentmem = -1;
    }

    /// Get a PDF member
    const PDFPtr member(int mem) const {
      // if (members.find(mem) == members.end())
      //   loadMember(mem);
      return members.find(mem)->second;
    }

    /// Get the currently active PDF member
    const PDFPtr activemember() const {
      return members.find(currentmem)->second;
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

  static const int PIDS[13] = { -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 };
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
    /// @todo Strip file extension for backward compatibility
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      ACTIVESETS[nset] = PDFSetHandler(setpath); //< @todo Will be wrong if a structured path is given
  }


  /// Load a PDF set by name
  void initpdfsetbynamem_(int& nset, const char* name, int namelength) {
    /// @todo Strip file extension for backward compatibility
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      ACTIVESETS[nset] = PDFSetHandler(name);
  }


  /// Load a PDF in current set.
  void initpdfm_(int& nset, int& nmember) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    ACTIVESETS[nset].loadMember(nmember);
  }


  /// Get xf(x) values for common partons from current PDF.
  void evolvepdfm_(int& nset, double& x, double& q , double* fxq) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // Evaluate for the LHAPDF5 standard partons
    for (size_t i = 0; i < 13; ++i) {
      try {
        fxq[i] = ACTIVESETS[nset].activemember()->xfxQ2(PIDS[i], x, q*q);
      } catch (const std::exception& e) {
        fxq[i] = 0;
      }
    }
  }


  /// Get photon PID xfx value from current PDF.
  void evolvepdfphotonm_(int& nset, double& x, double& q, double& fxq, double& photonparam) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    // Evaluate for the MRST2004QED photon flavor
    try {
      fxq = ACTIVESETS[nset].activemember()->xfxQ2(22, x, q*q);
    } catch (const std::exception& e) {
      fxq = 0;
    }
  }


  /// Set LHAPDF parameters
  void setlhaparm_(const char* par, int parlength) {
    // Do nothing for now
    /// @todo Can any Fortran LHAPDF params be usefully mapped?
  }


  /// @todo Add mapping of the xfx_ etc. functions for PYTHIA6 and friends. Use of LHAGLUE PDFLIB codes?


}
