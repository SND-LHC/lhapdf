from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "../../include/LHAPDF/Config.h" namespace "LHAPDF":
    cdef string version()

cdef extern from "../../include/LHAPDF/Paths.h" namespace "LHAPDF":
    cdef vector[string] _paths()
    cdef void setPaths(vector[string])
    cdef void pathsPrepend(string)
    cdef void pathsAppend(string)

cdef extern from "../../include/LHAPDF/PDF.h" namespace "LHAPDF":
    cdef cppclass PDF:
        double xfxQ(int, double, double) except +
        double xfxQ2(int, double, double) except +
        bool inRangeX(double) except +
        bool inRangeQ(double) except +
        bool inRangeQ2(double) except +
        bool inRangeXQ(double, double) except +
        bool inRangeXQ2(double, double) except +
        bool hasFlavor(int) except +
        vector[int] flavors()
        int numMembers() except +
        int memberID() except +
        int lhapdfID() except +
        string type() except +
        string setname() except +
        string description() except +
        string setdescription() except +
        string errorType() except +
        int qcdOrder() except +

cdef extern from "../../include/LHAPDF/Factories.h" namespace "LHAPDF":
    cdef PDF* mkPDF(string, int)
    cdef PDF* mkPDF(int)
