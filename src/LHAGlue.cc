#include "LHAPDF/Factories.h"
#include <sstream>
#include <iostream>

using namespace std;

namespace {
  static const int nPartons = 13;
  static int pid [nPartons] = { -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 };

  struct PDFSetHandler {
    PDFSetHandler(const string& name)
      : setname(setname)
    {    }

    /// Load a new PDF member
    void loadMember(int mem) {
      members[mem] = LHAPDF::loadPDF(setname, mem);
    }

    /// Actively delete a PDF member to save memory
    void unloadMember(int mem) {
      members.erase(mem);
    }

    /// Use smart pointers in the map to auto-delete
    typedef auto_ptr<PDF> PDFPtr;
    map<int, PDFPtr > members;
  };

  static map<string, PDFSetHandler> activesets;
}


extern "C" {

  /// LHAPDF version
  void getlhapdfversion_(char *s, size_t len) {
    LHAPDF::getVersion(len, s);
  }

  /// Load a PDF set
  void initpdfsetm_(int& nset, const char* setpath, int setpathlength) {
    /// @todo Strip file extension for backward compatibility
    if (wrapmap.find(nset) == wrapmap.end()) {
      /// @todo Adapt for lack of PDFSet layer... or reintroduce?
      wrapmap[nset].set = LHAPDF::PDFSet::load(setpath);
    }
  }

  /// Load a PDF set by name
  void initpdfsetbynamem_(int& nset, const char* name, int namelength) {
    /// @todo Strip file extension for backward compatibility
    if (wrapmap.find(nset) == wrapmap.end()) {
      /// @todo Actually map the name
      wrapmap[nset].set = LHAPDF::PDFSet::loadByName("MSTW2008lo90cl");
    }
  }

  /// Load a PDF in current set.
  void initpdfm_(int& nset, int& nmember) {
    std::map<int, PDFSetWrap>::iterator bundle = wrapmap.find(nset);
    if (bundle != wrapmap.end()) {
      bundle->second.pdf = &bundle->second.set->getMember(nmember);
    } else {
      throw UserError(nset);
    }
  }

  /// Get xf(x) values for common partons from current PDF.
  void evolvepdfm_(int& nset, double& x, double& q , double* fxq) {
    std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find(nset);
    // Check that set exists
    if ( bundle != wrapmap.end() ) {
      // Evaluate LHAPDF standard partons
      for (int i = 0; i < nPartons; i++) {
        try {
          fxq[i] = wrapmap[nset][nmem]->xfxQ2( pid[i], x, q*q );
        } catch ( std::exception e ) {
          fxq[i] = 0;
        }
      }
    }
  }

  /// Get photon PID xfx value from current PDF.
  void evolvepdfphotonm_(int& nset, double& x, double& q, double& fxq, double& photonparam) {
    std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );
    // Check that set exists
    if ( bundle != wrapmap.end() ) {
      fxq = bundle->second.pdf->xfxQ2( 22, x, q*q );
    } else {
      throw PDFSetNotLoadedException( nset );
    }
  }

  ///
  void setlhaparm_(const char* par, int parlength) {
    std::cout << "ERR: " << par << " SETLHAPARM" << std::endl;
  }

  /// @todo Add mapping of the xfx_ etc. functions for PYTHIA6 and friends

}
