from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "../../include/LHAPDF/Version.h" namespace "LHAPDF":
    cdef string version()

cdef extern from "../../include/LHAPDF/Paths.h" namespace "LHAPDF":
    cdef vector[string] paths()
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
        double alphasQ(double) except +
        double alphasQ2(double) except +
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
        void _print "print" () except + # TODO: add the second (verbosity) argument
        PDFSet& set()
        PDFInfo& info()
        #AlphaS& alphaS()

cdef extern from "../../include/LHAPDF/Info.h" namespace "LHAPDF":
    cdef cppclass Info:
        map[string,string] metadata()
        bool has_key(string)
        bool has_key_local(string)
        string get_entry(string) except +
        string get_entry(string, string) except +
        void set_entry(string, string)

cdef extern from "../../include/LHAPDF/Config.h" namespace "LHAPDF":
    cdef cppclass Config(Info.Info):
        pass
    cdef int verbosity()
    cdef void setVerbosity(int)

cdef extern from "../../include/LHAPDF/PDFSet.h" namespace "LHAPDF":
    cdef cppclass PDFSet(Info.Info):
        vector[PDF*] mkPDFs()
        PDF* mkPDF(int)
        size_t size() except +
        string name() except +
        string description() except +
        string errorType() except +
        void _print "print" () except + # TODO: map the second (verbosity) argument

cdef extern from "../../include/LHAPDF/PDFInfo.h" namespace "LHAPDF":
    cdef cppclass PDFInfo(Info.Info):
        pass #bool has_key(string)

cdef extern from "../../include/LHAPDF/Factories.h" namespace "LHAPDF":
    cdef Info& getConfig() #except +  # TODO: re-enable when Cython refs+exceptions has been bugfixed
    cdef PDFSet& getPDFSet(string) #except +  # TODO: re-enable when Cython refs+exceptions has been bugfixed
    cdef vector[PDF*] mkPDFs(string) except +
    cdef PDF* mkPDF(string, int) except +
    cdef PDF* mkPDF(int) except +
    #cdef AlphaS* mkAlphaS(string) except +
