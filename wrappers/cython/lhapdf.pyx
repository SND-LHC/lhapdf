cimport clhapdf as c
from libcpp.string cimport string
from libcpp.vector cimport vector
from itertools import izip


cdef class PDF:
    """
    """
    cdef c.PDF* _ptr
    cdef set_ptr(self, c.PDF* ptr):
        self._ptr = ptr

    def xfxQ2(self, pid, x, q2):
        try:
            return [self._ptr.xfxQ2(pid, x, q2) for x, q2 in izip(x, q2)]
        except TypeError:
            return self._ptr.xfxQ2(pid, x, q2)

    def xfxQ(self, pid, x, q):
        try:
            return [self._ptr.xfxQ(pid, x, q) for x, q in izip(x, q)]
        except TypeError:
            return self._ptr.xfxQ(pid, x, q)

    def inRangeQ(self, q):
        return self._ptr.inRangeQ(q)

    def inRangeQ2(self, q2):
        return self._ptr.inRangeQ2(q2)

    def inRangeX(self, x):
        return self._ptr.inRangeX(x)

    def inRangeXQ(self, x, q):
        return self._ptr.inRangeXQ(x, q)

    def inRangeXQ2(self, x, q2):
        return self._ptr.inRangeXQ2(x, q2)

    def hasFlavor(self, pid):
        return self._ptr.hasFlavor(pid)

    @property
    def memberid(self):
        return self._ptr.memberid()

    def __dealloc__(self):
        del self._ptr


def mkPDF(char* setname, int memid):
    cdef PDF obj = PDF.__new__(PDF)
    obj.set_ptr(c.mkPDF(string(setname), memid))
    return obj



# def do_something(char *name):
#     cdef c.PDFGrid *grid = c.PDFGrid_load(string(name))
#          print grid.xfxQ2(1, 0.3, 1000)

# cdef knotsToList(vector[double] knots):
#     outputKnots = []
#     cdef size_t i
#     for i in range(knots.size()):
#         outputKnots.append(knots[i])
#     return outputKnots


# cdef class PDFGrid(PDF):
#     """
#     """
#     cdef c.PDFGrid *_PDFGrid(self):
#         return <c.PDFGrid *> self._ptr

#     def inPhysicalRangeQ2(self, x, q2):
#         return self._PDFGrid().inPhysicalRangeQ2(x, q2)

#     def lookupClosestLow(self, x, q2, xid, q2id):
#         cdef c.Idx_t xi, q2i
#         xi, q2i = xid, q2id
#         self._PDFGrid().lookupClosestLow(x, q2, xi, q2i)
#         return xi, q2i

#     def lookupClosest(self, x, q2, xid, q2id):
#         cdef c.Idx_t xi, q2i
#         xi, q2i = xid, q2id
#         self._PDFGrid().lookupClosest(x, q2, xi, q2i)
#         return xi, q2i

#     #def isValidIndex(self, xid, q2id):
#     #    return self._PDFGrid().isValidIndex(xid, q2id)


#     @property
#     def xKnots(self):
#         return knotsToList(self._PDFGrid().getXKnots())

#     @property
#     def q2Knots(self):
#         return knotsToList(self._PDFGrid().getQ2Knots())

#     def setInterpolator(self, Interpolator interpolator):
#         self._PDFGrid().setInterpolator(interpolator._ptr)

#     def setDefaultInterpolator(self):
#         self._PDFGrid().setDefaultInterpolator()

#     def setExtrapolator(self, Extrapolator extrapolator):
#         self._PDFGrid().setExtrapolator(extrapolator._ptr)

#     def setDefaultExtrapolator(self):
#         self._PDFGrid().setDefaultExtrapolator()

#     def hasInterpolator(self):
#         return self._PDFGrid().hasInterpolator()

#     def hasExtrapolator(self):
#         return self._PDFGrid().hasExtrapolator()

#     @staticmethod
#     def load(char *fname):
#         cdef PDFGrid obj = PDFGrid.__new__(PDFGrid)
#         obj.set_ptr(c.PDFGrid_load(string(fname)))
#         return obj

#     @staticmethod
#     def loadByName(char *name):
#         cdef PDFGrid obj = PDFGrid.__new__(PDFGrid)
#         obj.set_ptr(c.PDFGrid_loadByName(string(name)))
#         return obj

#     def __dealloc__(self):
#         del self._ptr


# cdef class Interpolator:
#     """
#     """

#     cdef c.Interpolator *_ptr

#     cdef set_ptr(self, c.Interpolator *ptr):
#         self._ptr = ptr

#     def bind(self, PDFGrid pdfgrid):
#         self._ptr.bind(pdfgrid._PDFGrid()[0])

#     def unbind(self, PDFGrid pdfgrid):
#         self._ptr.unbind(pdfgrid._PDFGrid()[0])

#     def interpolateQ(self, PDFGrid pdfgrid, pid, x, q):
#         self._ptr.interpolateQ(pdfgrid._PDFGrid()[0], pid, x, q)

#     def interpolateQ2(self, PDFGrid pdfgrid, pid, x, q2):
#         self._ptr.interpolateQ2(pdfgrid._PDFGrid()[0], pid, x, q2)

#     def __dealloc__(self):
#         del self._ptr

#     #TODO: static methods


# cdef class Extrapolator:
#     """
#     """

#     cdef c.Extrapolator *_ptr

#     cdef set_ptr(self, c.Extrapolator *ptr):
#         self._ptr = ptr

#     def bind(self, PDFGrid pdfgrid):
#         self._ptr.bind(pdfgrid._PDFGrid()[0])

#     def unbind(self, PDFGrid pdfgrid):
#         self._ptr.unbind(pdfgrid._PDFGrid()[0])

#     def extrapolateQ(self, PDFGrid pdfgrid, pid, x, q):
#         self._ptr.extrapolateQ(pdfgrid._PDFGrid()[0], pid, x, q)

#     def extrapolateQ2(self, PDFGrid pdfgrid, pid, x, q2):
#         self._ptr.extrapolateQ2(pdfgrid._PDFGrid()[0], pid, x, q2)

#     def __dealloc__(self):
#         del self._ptr


# cdef class PDFSet:
#     """
#     """

#     cdef c.PDFSet *_ptr

#     cdef set_ptr(self, c.PDFSet *ptr):
#         self._ptr = ptr

#     #@property
#     #def member(self, member):
#     #    cdef PDF obj = PDF.__new__(PDF)
#     #    obj.set_ptr(self._ptr.getMember(member))
#     #    return obj

#     def alphaS(self, Q2):
#         return self._ptr.alphaS(Q2)

#     @staticmethod
#     def load(char *fname):
#         cdef PDFSet obj = PDFSet.__new__(PDFSet)
#         obj.set_ptr(c.PDFSet_load(string(fname)))
#         return obj

#     @staticmethod
#     def loadByName(char *name):
#         cdef PDFSet obj = PDFSet.__new__(PDFSet)
#         obj.set_ptr(c.PDFSet_loadByName(string(name)))
#         return obj


# def createInterpolator(char *interpolatorName):
#     cdef Interpolator obj = Interpolator.__new__(Interpolator)
#     obj.set_ptr(c.createInterpolator(string(interpolatorName)))
#     return obj

# def createExtrapolator(char *extrapolatorName):
#     cdef Extrapolator obj = Extrapolator.__new__(Extrapolator)
#     obj.set_ptr(c.createExtrapolator(string(extrapolatorName)))
#     return obj
