from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint32_t

ctypedef int PID_t
ctypedef double* PIDdata
ctypedef map[PID_t,PIDdata] PIDmap
ctypedef vector[double] AxisKnots
ctypedef double X_t
ctypedef double Q_t
ctypedef double Q2_t
ctypedef size_t Member_t
ctypedef size_t Idx_t

cdef extern from "../../include/PDF.h":
    cdef cppclass PDF:
        double xfxQ(PID_t, X_t, Q_t) except +
        double xfxQ2(PID_t, X_t, Q2_t) except +
        void xfxQ(PID_t, int, X_t, Q_t, double*) except +
        void xfxQ2(PID_t, int, X_t, Q2_t, double*) except + 
        bool inRangeQ(X_t, Q_t) except +
        bool inRangeQ2(X_t, Q2_t) except +
        bool hasPID(PID_t) except +
        Member_t getMemberNumber() except +

cdef extern from "../../include/PDFGrid.h":
    cdef cppclass PDFGrid(PDF):
        bool inPhysicalRangeQ2(X_t, Q2_t)
        void lookupClosestLow(X_t, Q2_t, Idx_t &, Idx_t &)
        void lookupClosest(X_t, Q2_t, Idx_t &, Idx_t &)
        bool isValidIndex(Idx_t, Idx_t)
        AxisKnots getXKnots()
        AxisKnots getQ2Knots()
        PIDdata getPIDData(PID_t)
        void setInterpolator(Interpolator*)
        void setDefaultInterpolator()
        void setExtrapolator(Extrapolator*)
        void setDefaultExtrapolator()
        bool hasInterpolator()
        bool hasExtrapolator()

    cdef PDFGrid *PDFGrid_load "PDFGrid::load" (string)   
    cdef PDFGrid *PDFGrid_loadByName "PDFGrid::loadByName" (string)

cdef extern from "../../include/Interpolator.h":
    cdef cppclass Interpolator:
        void bind(PDFGrid)
        void unbind(PDFGrid)
        double interpolateQ(PDFGrid, PID_t, X_t, Q_t)
        double interpolateQ2(PDFGrid, PID_t, X_t, Q2_t) 
        void interpolateQ(PDFGrid, PID_t, uint32_t, X_t, Q_t, double)  
        void interpolateQ2(PDFGrid, PID_t, uint32_t, X_t, Q_t, double)
    
    cdef void Interpolator_bind "Interpolator::bind" (PDFGrid)
    cdef void Interpolator_unbind "Interpolator::unbind" (PDFGrid)
    cdef double Interpolator_interpolateQ "Interpolator::interpolateQ" (
        PDFGrid, PID_t, X_t, Q_t)

# cdef void Interpolator_interpolateQ "Interpolator::interpolateQ" (PDFGrid,
#                                                                       PID_t,
#                                                                       uint32_t,
#                                                                       X_t, Q_t,
#                                                                       double)
# cdef void Interpolator_interpolateQ2 "Interpolator::interpolateQ2" (PDFGrid,
#                                                                         PID_t,
#                                                                         uint32_t,
#                                                                         X_t,
#                                                                         Q2_t,
#                                                                        double)

cdef extern from "../../include/Extrapolator.h":
    cdef cppclass Extrapolator:
        void bind(PDFGrid)
        void unbind(PDFGrid)
        double extrapolateQ(PDFGrid, PID_t, X_t, Q_t)
        double extrapolateQ2(PDFGrid, PID_t, X_t, Q2_t) 
        void extrapolateQ(PDFGrid, PID_t, uint32_t, X_t, Q_t, double)  
        void extrapolateQ2(PDFGrid, PID_t, uint32_t, X_t, Q_t, double)

cdef extern from "../../include/PDFSet.h":
    cdef cppclass PDFSet:
        PDF *getMember(Member_t)
        double alphaS(Q2_t)

    cdef PDFSet *PDFSet_load "PDFSet::load" (string)
    cdef PDFSet *PDFSet_loadByName "PDFSet::loadByName" (string)

cdef extern from "../../include/bicubic.h":
    cdef cppclass BicubicInterpolator(Interpolator):
        pass

cdef extern from "../../include/bilinear.h":
    cdef cppclass BilinearInterpolator(Interpolator):
        pass

cdef extern from "../../include/factory.h":
    cdef Interpolator* createInterpolator(string)
    cdef Extrapolator* createExtrapolator(string)

cdef extern from "../../include/nearest.h":
    cdef cppclass NearestExtrapolator(Extrapolator):
        pass

cdef extern from "../../include/version.h":
    cdef int getVersion(int, string)
