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
    def qcdOrder(self):
        "Max number of loops involved in this PDF's evolution."
        return self._ptr.qcdOrder()

    @property
    def xMin(self):
        "Minimum valid value of x to be used with this PDF"
        return self._ptr.xMin()

    @property
    def xMax(self):
        "Maximum valid value of x to be used with this PDF"
        return self._ptr.xMax()

    @property
    def q2Min(self):
        "Minimum valid value of x to be used with this PDF"
        return self._ptr.q2Min()

    @property
    def q2Max(self):
        "Maximum valid value of x to be used with this PDF"
        return self._ptr.q2Max()

    def alphasQ(self, q):
        "Return alpha_s at q"
        return self._ptr.alphasQ(q)

    def alphasQ2(self, q2):
        "Return alpha_s at q2"
        return self._ptr.alphasQ2(q2)

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

    cdef _set(self):
        cdef PDFSet obj = PDFSet.__new__(PDFSet)
        obj.set_ptr(&self._ptr.set())
        return obj

    def set(self):
        "Return the corresponding PDFSet"
        return self._set()

    cdef _info(self):
        cdef PDFInfo obj = PDFInfo.__new__(PDFInfo)
        obj.set_ptr(&self._ptr.info())
        return obj

    def info(self):
        "Return the corresponding PDFInfo"
        return self._info()

    def _print(self):
        "Print a short summary to stdout"
        self._ptr._print()


cdef class Info:
    """\
    Class that handles the parsing of PDF set metadata in the .info file.
    """
    cdef c.Info* _ptr
    cdef set_ptr(self, c.Info* ptr):
        self._ptr = ptr

    def metadata(self):
        "Return the metadata in the .info file"
        return self._ptr.metadata()

    def has_key(self, key):
        "Return whether or not metadata for this key exists"
        return self._ptr.has_key(key)

    def has_key_local(self, key):
        "Returns whether or not metadata for this key exists at a local level (config/set/member)"
        return self._ptr.has_key_local(key)

    # def get_entry(self, key):
    #     "Returns metadata entry for this key"
    #     return self._ptr.get_entry(key)

    def get_entry(self, key, fallback=None):
        "Returns metadata entry for this key if it exists, otherwise returns a fallback value"
        rtn = self._ptr.get_entry(key, str(fallback))
        return rtn if str(rtn) != str(fallback) else fallback

    def set_entry(self, key, value):
        "Set a metadata key"
        self._ptr.set_entry(key, str(value))


cdef class PDFSet:
    """\
    A collection of PDFs with related fits, most typically a central PDF and a
    set of extra ones representing different aspects of systematic errors in the
    fit.
    """
    cdef c.PDFSet* _ptr
    cdef set_ptr(self, c.PDFSet* ptr):
        self._ptr = ptr

    def __dealloc__(self):
        pass

    def __len__(self):
        "The total number of members in this set."
        return self._ptr.size()

    @property
    def size(self):
        "The total number of members in this set."
        return self._ptr.size()

    @property
    def name(self):
        "Name of this PDF's containing set."
        return self._ptr.name()

    @property
    def description(self):
        "Description of this PDF's set."
        return self._ptr.description()

    @property
    def errorType(self):
        "Type of error treatment in this PDF's set."
        return self._ptr.errorType()

    def mkPDF(self, mem):
        cdef c.PDF* ptr = self._ptr.mkPDF(mem)
        cdef PDF obj
        obj = PDF.__new__(PDF)
        obj.set_ptr(ptr)
        return obj

    def mkPDFs(self):
        cdef vector[c.PDF*] ptrs = self._ptr.mkPDFs()
        cdef PDF obj
        objs = []
        for ptr in ptrs:
            obj = PDF.__new__(PDF)
            obj.set_ptr(ptr)
            objs.append(obj)
        return objs

    def metadata(self):
        "Return the metadata in the .info file"
        return self._ptr.metadata()

    def has_key(self, key):
        "Return whether or not metadata for this key exists"
        return self._ptr.has_key(key)

    def has_key_local(self, key):
        "Returns whether or not metadata for this key exists at a local level (config/set/member)"
        return self._ptr.has_key_local(key)

    # def get_entry(self, key):
    #     "Returns metadata entry for this key"
    #     return self._ptr.get_entry(key)

    def get_entry(self, key, fallback=None):
        "Returns metadata entry for this key if it exists, otherwise returns a fallback value"
        rtn = self._ptr.get_entry(key, str(fallback))
        return rtn if str(rtn) != str(fallback) else fallback

    def _print(self):
        "Print a short summary to stdout"
        self._ptr._print()



cdef class PDFInfo:
    """\
    A class handling the metadata that defines a given PDF.
    """

    cdef c.PDFInfo* _ptr
    cdef set_ptr(self, c.PDFInfo* ptr):
        self._ptr = ptr

    def metadata(self):
        "Return the metadata in the .info file"
        return self._ptr.metadata()

    def has_key(self, key):
        "Return whether or not metadata for this key exists"
        return self._ptr.has_key(key)

    def has_key_local(self, key):
        "Returns whether or not metadata for this key exists at a local level (config/set/member)"
        return self._ptr.has_key_local(key)

    # def get_entry(self, key):
    #     "Returns metadata entry for this key"
    #     return self._ptr.get_entry(key)

    def get_entry(self, key, fallback=None):
        "Returns metadata entry for this key if it exists, otherwise returns a fallback value"
        rtn = self._ptr.get_entry(key, str(fallback))
        return rtn if str(rtn) != str(fallback) else fallback



def getConfig():
    """Factory function to get the global config object."""
    cdef c.Info* ptr = &c.getConfig()
    cdef Info obj = Info.__new__(Info)
    obj.set_ptr(ptr)
    return obj


def getPDFSet(setname):
    """Factory function to get the specified PDF set."""
    cdef c.PDFSet* ptr = &c.getPDFSet(setname)
    cdef PDFSet obj = PDFSet.__new__(PDFSet)
    obj.set_ptr(ptr)
    return obj

def mkPDFs(setname):
    """Factory function to make all the PDF objects in the specified set."""
    cdef vector[c.PDF*] ptrs = c.mkPDFs(setname)
    cdef PDF obj
    objs = []
    for ptr in ptrs:
        obj = PDF.__new__(PDF)
        obj.set_ptr(ptr)
        objs.append(obj)
    return objs


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



# TODO: map AlphaS



def version():
    "Return the LHAPDF library version."
    return c.version()

__version__ = version()



def verbosity():
    "Get the main verbosity level of the LHAPDF system: 0 = quiet, 2 = loud"
    return c.verbosity()

def setVerbosity(vlevel):
    "Set the main verbosity level of the LHAPDF system: 0 = quiet, 2 = loud"
    c.setVerbosity(vlevel)


def availablePDFSets():
    "Get the names of all the available PDF sets on this system."
    return c.availablePDFSets()


def paths():
    "Return the list of current PDF data search paths."
    return c._paths()

def setPaths(newpaths):
    "Set the list of current PDF data search paths."
    c.setPaths(newpaths)

def pathsPrepend(newpath):
    "Prepend to the list of current PDF data search paths."
    c.pathsPrepend(newpath)

def pathsAppend(newpath):
    "Append to the list of current PDF data search paths."
    c.pathsAppend(newpath)
