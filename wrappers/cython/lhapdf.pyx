cimport clhapdf as c
from libcpp.string cimport string

# def do_something(char *name):
#     cdef c.PDFGrid *grid = c.PDFGrid_load(string(name))

#          print grid.xfxQ2(1, 0.3, 1000)

cdef class PDF:
    """
    """
    cdef c.PDF *_ptr
    
    def xfxQ2(self, PID_t, X_t, Q2_t):
	return self._ptr.xfxQ2(PID_t, X_t, Q2_t)

    def xfxQ(self, PID_t, X_t, Q_t):
	return self._ptr.xfxQ(PID_t, X_t, Q_t)

    def __dealloc__(self):
	del self._ptr 

cdef class PDFGrid:
    """
    """
    cdef c.PDFGrid *_ptr

    def *PDFGrid_load(string):
	self._ptr.load(string)
