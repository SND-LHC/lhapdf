from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "../../include/LHAPDF/Config.h" namespace "LHAPDF":
    cdef int getVersion(int, string)

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
        # TODO: vector<int> flavors()
        int memberid() except +

cdef extern from "../../include/LHAPDF/Factories.h" namespace "LHAPDF":
    cdef PDF* mkPDF(string, int)
    cdef PDF* mkPDF(int)


# cdef extern from "../../include/LHAPDF/GridPDF.h":
#     cdef cppclass GridPDF(PDF):
#         bool inPhysicalRangeQ2(double, double)
#         void lookupClosestLow(double, double, Idx_t &, Idx_t &)
#         void lookupClosest(double, double, Idx_t &, Idx_t &)
#         AxisKnots getXKnots()
#         AxisKnots getQ2Knots()
#         PIDdata getPIDData(int)
#         void setInterpolator(Interpolator*)
#         void setDefaultInterpolator()
#         void setExtrapolator(Extrapolator*)
#         void setDefaultExtrapolator()
#         bool hasInterpolator()
#         bool hasExtrapolator()

#     cdef PDFGrid *PDFGrid_load "PDFGrid::load" (string)
#     cdef PDFGrid *PDFGrid_loadByName "PDFGrid::loadByName" (string)

# cdef extern from "../../include/Interpolator.h":
#     cdef cppclass Interpolator:
#         void bind(PDFGrid)
#         void unbind(PDFGrid)
#         double interpolateQ(PDFGrid, int, double, double)
#         double interpolateQ2(PDFGrid, int, double, double)
#         void interpolateQ(PDFGrid, int, uint32_t, double, double, double)
#         void interpolateQ2(PDFGrid, int, uint32_t, double, double, double)

#     cdef void Interpolator_bind "Interpolator::bind" (PDFGrid)
#     cdef void Interpolator_unbind "Interpolator::unbind" (PDFGrid)
#     cdef double Interpolator_interpolateQ "Interpolator::interpolateQ" (
#         PDFGrid, int, double, double)

# cdef void Interpolator_interpolateQ "Interpolator::interpolateQ" (PDFGrid,
#                                                                       int,
#                                                                       uint32_t,
#                                                                       double, double,
#                                                                       double)
# cdef void Interpolator_interpolateQ2 "Interpolator::interpolateQ2" (PDFGrid,
#                                                                         int,
#                                                                         uint32_t,
#                                                                         double,
#                                                                         double,
#                                                                        double)

# cdef extern from "../../include/Extrapolator.h":
#     cdef cppclass Extrapolator:
#         void bind(PDFGrid)
#         void unbind(PDFGrid)
#         double extrapolateQ(PDFGrid, int, double, double)
#         double extrapolateQ2(PDFGrid, int, double, double)
#         void extrapolateQ(PDFGrid, int, uint32_t, double, double, double)
#         void extrapolateQ2(PDFGrid, int, uint32_t, double, double, double)


# cdef extern from "../../include/Bicubic.h":
#     cdef cppclass BicubicInterpolator(Interpolator):
#         pass

# cdef extern from "../../include/Bilinear.h":
#     cdef cppclass BilinearInterpolator(Interpolator):
#         pass

# cdef extern from "../../include/factory.h":
#     cdef Interpolator* createInterpolator(string)
#     cdef Extrapolator* createExtrapolator(string)

# cdef extern from "../../include/Nearest.h":
#     cdef cppclass NearestExtrapolator(Extrapolator):
#         pass

# cdef extern from "../../include/version.h":
#     cdef int getVersion(int, string)
