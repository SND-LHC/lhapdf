#! /usr/bin/env python

## Python LHAPDF6 usage example for PDF uncertainty code (G. Watt, 24/04/2014)

import lhapdf
import numpy as np

x = 0.1
q = 100.0
pset = lhapdf.getPDFSet("CT10nnlo")
#pset = lhapdf.getPDFSet("NNPDF23_nnlo_as_0118")
pdfs = pset.mkPDFs()
nmem = pset.size - 1

print
print "Error type = ", pset.errorType
print "Error conf level = ", pset.errorConfLevel

## Fill vectors xgAll and xuAll using all PDF members.
xgAll = np.empty([pset.size])
xuAll = np.empty([pset.size])
for imem in xrange(pset.size):
    xgAll[imem] = pdfs[imem].xfxQ(21, x, q)
    xuAll[imem] = pdfs[imem].xfxQ(2, x, q)
    # print imem, xgAll[imem], xuAll[imem]

## Calculate 1-sigma PDF uncertainty on gluon distribution.
unc = pset.uncertainty(xgAll)
print "xg = %.3g + %.3g - %.3g (+- %.3g)" % (unc.central, unc.errplus, unc.errminus, unc.errsymm)
print "Scale = ", unc.scale

## Calculate 1-sigma PDF uncertainty on up-quark distribution.
unc = pset.uncertainty(xuAll)
print "xu = %.3g + %.3g - %.3g (+- %.3g)" % (unc.central, unc.errplus, unc.errminus, unc.errsymm)
print "Scale = ", unc.scale

## Calculate PDF correlation between gluon and up-quark.
corr = pset.correlation(xgAll, xuAll)
print "Correlation = ",corr

## Generate random values from Hessian best-fit and eigenvector values.
## Random numbers are just set to zero here for testing purposes.
#randoms = np.empty([nmem/2])
randoms = [0.0 for ir in xrange(nmem/2)]
xgrand = pset.randomValueFromHessian(xgAll, randoms)
xurand = pset.randomValueFromHessian(xuAll, randoms)
print "Random: xg = %.3g, xu = %.3g" % (xgrand, xurand)
