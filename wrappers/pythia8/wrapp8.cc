#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/Exceptions.h"
#include <sstream>
#include <exception>
#include <iostream>

using namespace LHAPDF;


/// PDFSetNotLoadedException is thrown from the Pythia8 wrapper when the selected PDFSet is not loaded.
class PDFSetNotLoadedException : public Exception {
public:
  PDFSetNotLoadedException( const int Set )
    : set(Set) {}

  const char* what() const throw() {
    std::stringstream message;
    message << "Trying to access set: " << set << " which is not loaded.";
    return message.str().c_str();
  }

private:
  const int set;
};


/// A structure that is required to set a current PDF member similar to the original
/// LHAPDF interface.
struct PDFSetWrap {
  PDFSet* set;
  PDF* pdf;
};

static const int nPartons = 13;
static int pid [nPartons] = {-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6};

static std::map<int, PDFSetWrap> wrapmap;

extern "C" {
  /// LHAPDF Version
  void getlhapdfversion_(char *s, size_t len) {
    getVersion(len, s);
  }

  /// Load a PDFSet.
  void initpdfsetm_(int& nset, const char* setpath, int setpathlength) {
    if (wrapmap.find(nset) == wrapmap.end()) {
      //Set does not exist
      wrapmap[nset].set = PDFSet::load( setpath );
    }
  }

  /// Load a PDFSet by name.
  void initpdfsetbynamem_(int& nset, const char* name, int namelength) {
    if ( wrapmap.find( nset ) == wrapmap.end() ) {
      //Set does not exist
      //wrapmap[nset].set = PDFSet::loadByName( name );
      wrapmap[nset].set = PDFSet::loadByName( "MSTW2008lo90cl" );

    }
  }

  /// Load a PDF in current PDFSet.
  void initpdfm_(int& nset, int& nmember) {
    std::map<int, PDFSetWrap>::iterator bundle = wrapmap.find( nset );

    //Check if PDFSet has been loaded
    if ( bundle != wrapmap.end() ) {
      bundle->second.pdf = &bundle->second.set->getMember( nmember );
    }
    else {
      throw PDFSetNotLoadedException( nset );
    }
  }

  /// Get xfx values for common partons from current PDF.
  void evolvepdfm_(int& nset, double& x, double& q , double* fxq) {
    std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );

    //Check that set exists
    if ( bundle != wrapmap.end() ) {
      //Evaluate LHAPDF standard partons
      for (int i = 0; i < nPartons; i++) {
        try {
          fxq[i] = bundle->second.pdf->xfxQ2( pid[i], x, q*q );
        } catch ( std::exception e ) {
          fxq[i] = 0;
        }
      }
    }
  }

  /// Get photon PID xfx value from current PDF.
  void evolvepdfphotonm_(int& nset, double& x, double& q, double& fxq, double& photonparam) {
    std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );

    //Check that set exists
    if ( bundle != wrapmap.end() ) {
      fxq = bundle->second.pdf->xfxQ2( 22, x, q*q );
    }
    else {
      throw PDFSetNotLoadedException( nset );
    }
  }

  ///
  void setlhaparm_(const char* par, int parlength) {
    std::cout << "ERR: " << par << " SETLHAPARM" << std::endl;
  }
}
