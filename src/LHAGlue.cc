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
    // typedef auto_ptr<LHAPDF::PDF> PDFPtr;
    typedef LHAPDF::PDF* PDFPtr;

    PDFSetHandler() { } //< It is stored in a map so we need one of these...

    PDFSetHandler(const string& name)
      : currentmem(-1), setname(setname)
    {    }

    PDFSetHandler(int lhaid) {
      pair<string,int> set_mem = LHAPDF::lookupPDF(lhaid);
      setname = set_mem.first;
      loadMember(set_mem.second);
    }

    /// Load a new PDF member
    void loadMember(int mem) {
      /// @todo What if it's already loaded? Don't leak memory
      members[mem] = LHAPDF::mkPDF(setname, mem);
      currentmem = mem;
    }

    /// Actively delete a PDF member to save memory
    void unloadMember(int mem) {
      members.erase(mem);
      currentmem = -1;
    }

    /// Get a PDF member
    const PDFPtr member(int mem) {
      return members[mem];
    }

    /// Get the currently active PDF member
    const PDFPtr activemember() {
      return members[currentmem];
    }

    /// The currently active member in this set
    int currentmem;

    /// Name of this set
    string setname;

    /// Map of pointers to selected member PDFs
    map<int, PDFPtr > members;
  };

  static const int PIDS[13] = { -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 };
  static map<int, PDFSetHandler> ACTIVESETS;
  // static int ACTIVESET, ACTIVEMEM;

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
    // ACTIVESET = nset
    // ACTIVEMEM = 0
  }


  /// Load a PDF set by name
  void initpdfsetbynamem_(int& nset, const char* name, int namelength) {
    /// @todo Strip file extension for backward compatibility
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      ACTIVESETS[nset] = PDFSetHandler(name);
    // ACTIVESET = nset
    // ACTIVEMEM = 0
  }


  /// Load a PDF in current set.
  void initpdfm_(int& nset, int& nmember) {
    if (ACTIVESETS.find(nset) == ACTIVESETS.end())
      throw LHAPDF::UserError("Trying to use LHAGLUE set #" + LHAPDF::to_str(nset) + " but it is not initialised");
    ACTIVESETS[nset].loadMember(nmember);
    // ACTIVESET = nset;
    // ACTIVEMEM = nmember;
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


  /// @todo Add mapping of the xfx_ etc. functions for PYTHIA6 and friends


}
