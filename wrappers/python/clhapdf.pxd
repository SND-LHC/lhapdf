from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "../../include/LHAPDF/Version.h" namespace "LHAPDF":
    cdef string version()

cdef extern from "../../include/LHAPDF/Paths.h" namespace "LHAPDF":
    cdef vector[string] _paths()
    cdef void setPaths(vector[string])
    cdef void pathsPrepend(string)
    cdef void pathsAppend(string)
    cdef vector[string] availablePDFSets()

cdef extern from "../../include/LHAPDF/PDF.h" namespace "LHAPDF":
    cdef cppclass PDF:
        double xfxQ(int, double, double) except +
        double xfxQ2(int, double, double) except +
        map[int,double] xfxQ(double, double) except +
        map[int,double] xfxQ2(double, double) except +
        double xMin()
        double xMax()
        double q2Min()
        double q2Max()
        bool inRangeX(double) except +
        bool inRangeQ(double) except +
        bool inRangeQ2(double) except +
        bool inRangeXQ(double, double) except +
        bool inRangeXQ2(double, double) except +
        bool hasFlavor(int) except +
        vector[int] flavors()
        int memberID() except +
        int lhapdfID() except +
        string description() except +
        string type() except +
        int qcdOrder() except +
        #PDFSet& set()
        #PDFInfo& info()
        #AlphaS& alphaS()

cdef extern from "../../include/LHAPDF/PDFSet.h" namespace "LHAPDF":
    cdef cppclass PDFSet:
        vector[PDF*] mkPDFs()
        PDF* mkPDF(int)
        # int numMembers() except +
        size_t size() except +
        string name() except +
        string description() except +
        string errorType() except +

cdef extern from "../../include/LHAPDF/Factories.h" namespace "LHAPDF":
    #cdef Info& getConfig()
    cdef PDFSet& getPDFSet(string)
    cdef vector[PDF*] mkPDFs(string)
    cdef PDF* mkPDF(string, int)
    cdef PDF* mkPDF(int)
    #cdef AlphaS* mkAlphaS(string)
