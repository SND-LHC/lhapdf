#cython: embedsignature=True

cimport clhapdf as c
from libcpp.string cimport string
from libcpp.vector cimport vector
from itertools import izip


cdef class PDF:
    """\
    A parton density function for in general several parton flavours,
    i.e. one member of a PDF set.
    """
    cdef c.PDF* _ptr
    cdef set_ptr(self, c.PDF* ptr):
        self._ptr = ptr

    def __dealloc__(self):
        del self._ptr

    @property
    def memberID(self):
        "The PDF set member number of this PDF."
        return self._ptr.memberID()

    @property
    def lhapdfID(self):
        "The LHAPDF ID number of this PDF member."
        return self._ptr.lhapdfID()

    @property
    def numMembers(self):
        "The total number of members in this set."
        return self._ptr.numMembers()

    # TODO: Need another name than "type" in Python?
    @property
    def type(self):
        "The type of PDF member, e.g. central, error."
        return self._ptr.type()

    @property
    def description(self):
        "Description of this PDF member."
        return self._ptr.description()

    @property
    def setdescription(self):
        "Description of this PDF's set."
        return self._ptr.setdescription()

    @property
    def setname(self):
        "Name of this PDF's containing set."
        return self._ptr.setname()

    @property
    def errorType(self):
        "Type of error treatment in this PDF's set."
        return self._ptr.errorType()

    @property
    def qcdOrder(self):
        "Max number of loops involved in this PDF's evolution."
        return self._ptr.qcdOrder()


    # TODO: alphaS


    def xfxQ(self, pid, x, q):
        # TODO: allow 2-arg version without PID which returns a dict for all flavours
        "Return the PDF xf(x,Q) value for the given parton ID, x, and Q."
        try:
            return [self._ptr.xfxQ(pid, x, q) for x, q in izip(x, q)]
        except TypeError:
            return self._ptr.xfxQ(pid, x, q)

    def xfxQ2(self, pid, x, q2):
        # TODO: allow 2-arg version without PID which returns a dict for all flavours
        "Return the PDF xf(x,Q2) value for the given parton ID, x, and Q2."
        try:
            return [self._ptr.xfxQ2(pid, x, q2) for x, q2 in izip(x, q2)]
        except TypeError:
            return self._ptr.xfxQ2(pid, x, q2)

    def inRangeQ(self, q):
        "Check if the specified Q value is in the unextrapolated range of this PDF."
        return self._ptr.inRangeQ(q)

    def inRangeQ2(self, q2):
        "Check if the specified Q2 value is in the unextrapolated range of this PDF."
        return self._ptr.inRangeQ2(q2)

    def inRangeX(self, x):
        "Check if the specified x value is in the unextrapolated range of this PDF."
        return self._ptr.inRangeX(x)

    def inRangeXQ(self, x, q):
        "Check if the specified x and Q values are in the unextrapolated range of this PDF."
        return self._ptr.inRangeXQ(x, q)

    def inRangeXQ2(self, x, q2):
        "Check if the specified x and Q2 values are in the unextrapolated range of this PDF."
        return self._ptr.inRangeXQ2(x, q2)

    def flavors(self):
        "Return the list of parton IDs supported by this PDF."
        # TODO: Use Cython >= 0.17 STL type coercion when available
        cdef vector[int] flavs = self._ptr.flavors()
        return [flavs[i] for i in xrange(flavs.size())]

    def hasFlavor(self, pid):
        "Check if the specified parton ID is contained in this PDF."
        return self._ptr.hasFlavor(pid)

    # TODO: Map the rest of the metadata functions (including the generic metadata() -> str)


cdef mkPDF_setmem(char* setname, int memid):
    "Factory function to make a PDF object from the set name and member number."
    cdef PDF obj = PDF.__new__(PDF)
    obj.set_ptr(c.mkPDF(string(setname), memid))
    return obj

cdef mkPDF_lhaid(int lhaid):
    "Factory function to make a PDF object from the LHAPDF ID number."
    cdef PDF obj = PDF.__new__(PDF)
    obj.set_ptr(c.mkPDF(lhaid))
    return obj

def mkPDF(*args):
    """Factory function to make a PDF object from the set name and member number
    (2 args), or just the unique LHAPDF ID number for that member (1 arg)."""
    if len(args) == 1 and type(args[0]) == int:
        return mkPDF_lhaid(args[0])
    elif len(args) == 2 and type(args[0]) == str and type(args[1]) == int:
        return mkPDF_setmem(args[0], args[1])
    else:
        raise Exception("Unknown call signature")


def version():
    "Return the LHAPDF library version."
    return str(c.version().c_str())
__version__ = version()


def paths():
    "Return the list of current PDF data search paths."
    # TODO: Use Cython >= 0.17 STL type coercion when available
    cdef vector[string] p = c._paths()
    return [p[i].c_str() for i in xrange(p.size())]

cdef _setPaths(newpaths):
    # TODO: Use Cython >= 0.17 STL type coercion when available
    cdef vector[string] vect
    cdef char* cs
    for s in newpaths:
        cs = s
        vect.push_back(string(cs))
    c.setPaths(vect)
def setPaths(newpaths):
    "Set the list of current PDF data search paths."
    _setPaths(newpaths)

cdef _pathsPrepend(char* newpath):
    # TODO: Use Cython >= 0.17 STL type coercion when available
    cdef string s = string(newpath)
    c.pathsPrepend(s)
def pathsPrepend(newpath):
    "Prepend to the list of current PDF data search paths."
    _pathsPrepend(newpath)

cdef _pathsAppend(char* newpath):
    # TODO: Use Cython >= 0.17 STL type coercion when available
    cdef string s = string(newpath)
    c.pathsAppend(s)
def pathsAppend(newpath):
    "Append to the list of current PDF data search paths."
    _pathsAppend(newpath)
