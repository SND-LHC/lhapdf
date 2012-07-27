from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector

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
        double xfxQ2(PID_t, X_t, Q2_t)
        double xfxQ(PID_t, X_t, Q_t)

cdef extern from "../../include/PDFGrid.h":
    cdef cppclass PDFGrid(PDF):
        pass
    
    cdef PDFGrid *PDFGrid_load "PDFGrid::load" (string)   
