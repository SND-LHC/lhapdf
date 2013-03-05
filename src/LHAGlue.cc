#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"

using namespace std;

namespace { //< Unnamed namespace to restrict visibility to this file


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

    /// @brief Internal storage is a smart pointer to ensure deletion of created PDFs
    ///
    /// NB. std::auto_ptr cannot be stored in STL containers, hence we use
    /// boost::shared_ptr. std::unique_ptr is the nature replacement when C++11
    /// is globally available.
    typedef boost::shared_ptr<LHAPDF::PDF> PDFPtr;

    /// Default constructor
    PDFSetHandler() : currentmember(0)
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
      const int nextmem = (members.size() != 0) ? members.begin()->first : 0;
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
