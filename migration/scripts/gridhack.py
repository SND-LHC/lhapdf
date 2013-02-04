from ctypes import *

try:
    _lha = CDLL('libLHAPDF.so')
except:
    _lha = CDLL('libLHAPDF.dylib')

def get_grid():
    nxgrid = c_uint()
    nqgrid = c_uint()
    gridx = (c_double * 500)()
    gridq = (c_double * 500)()
    _lha.getgrid_(byref(nxgrid), byref(nqgrid), gridx, gridq)
    return gridx[:nxgrid.value], gridq[:nqgrid.value]
