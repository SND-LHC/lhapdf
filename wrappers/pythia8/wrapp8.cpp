#include <sstream>
#include <exception>

#include "version.h"
#include "PDFSet.h"
#include "PDFGrid.h"

/// PDFSetNotLoadedException is thrown from the Pythia8 wrapper when the selected PDFSet is not loaded.
class PDFSetNotLoadedException: public std::exception {
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

/// A structure that is required to set a current pdf member similar to the original
/// LHAPDF interface.
struct PDFSetWrap {
    PDFSet* set;
    PDFGrid* pdf;
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
            //wrapmap[nset].set = PDFSet::load( setpath );
            
            wrapmap[nset].set = 0;
        }
        else {
            //TODO: What do we do when it is loaded?
        }
    }
    
    /// Load a PDFSet by name.
    void initpdfsetbynamem_(int& nset, const char* name, int namelength) {	
        if ( wrapmap.find( nset ) == wrapmap.end() ) {
            //Set does not exist
            //wrapmap[nset].set = PDFSet::loadByName( name );
            
            wrapmap[nset].set = 0;
        }
        else {
            //TODO: What do we do when it is loaded?
        }
    }

    /// Load a PDF in current PDFSet.
    void initpdfm_(int& nset, int& nmember) {
        std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );
        
        //Check if PDFSet has been loaded
        if ( bundle != wrapmap.end() ) {
            //bundle->second.set->getMember( nmember );
            
            bundle->second.pdf = PDFGrid::load( "/Home/s0821167/Summer Project/Interpolator/lhapdf/Tests/MSTW2008/MSTW2008lo90cl.LHm" );
        }
        else {
            throw PDFSetNotLoadedException( int );
        }
    }

    /// Get xfx values for common partons from current PDF.
    void evolvepdfm_(int& nset, double& x, double& q , double* fxq) {
        std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );
        
        //Check that set exists
        if ( bundle != wrapmap.end() ) {
            //Evaluate LHAPDF standard partons
            for (int i = 0; i < nPartons; i++) {
                fxq[i] = bundle->pdf->xfxQ2( pid[i], x, q*q );
            }
        }
    }

    /// Get photon PID xfx value from current PDF.
    void evolvepdfphotonm_(int& nset, double& x, double& q, double& fxq, double& photonparam) {        
        std::map<int, PDFSetWrap>::const_iterator bundle = wrapmap.find( nset );
        
        //Check that set exists
        if ( bundle != wrapmap.end() ) {
            fxq = bundle->xfxQ2( 22, x, q*q );
        }
        else {
            throw PDFSetNotLoadedException();
        }
    }

    /// 
    void setlhaparm_(const char* par, int parlength) {
        std::cout << "ERR: " << par << " SETLHAPARM" << std::endl;
    }
}
