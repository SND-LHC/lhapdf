#! /usr/bin/env python

## Python LHAPDF6 usage example for PDF uncertainty code (G. Watt, 24/04/2014).
## Extended in July 2015 for ErrorType values ending in "+as".

import lhapdf

x = 0.1
q = 100.0
pset = lhapdf.getPDFSet("CT10nnlo") # ErrorType: hessian
#pset = lhapdf.getPDFSet("abm12lhc_5_nnlo") # ErrorType: symmhessian
#pset = lhapdf.getPDFSet("NNPDF30_nnlo_as_0118") # ErrorType: replicas
#pset = lhapdf.getPDFSet("NNPDF30_nnlo_nf_5_pdfas") # ErrorType: replicas+as
pdfs = pset.mkPDFs()
nmem = pset.size - 1

print
print "Error type = ", pset.errorType
print "Error conf level = ", pset.errorConfLevel
print

## Fill vectors xgAll and xuAll using all PDF members.
xgAll = [0.0 for i in xrange(pset.size)]
xuAll = [0.0 for i in xrange(pset.size)]
for imem in xrange(pset.size):
    xgAll[imem] = pdfs[imem].xfxQ(21, x, q)
    xuAll[imem] = pdfs[imem].xfxQ(2, x, q)
    #print imem, xgAll[imem], xuAll[imem]

## Calculate 1-sigma PDF uncertainty on gluon distribution.
unc = pset.uncertainty(xgAll)
print "xg = %.4g + %.4g - %.4g (+- %.4g)" % (unc.central, unc.errplus, unc.errminus, unc.errsymm)
print "Scale = %.4g" % unc.scale
if pset.errorType.endswith("+as"):
    ## If the ErrorType ends in "+as" the last two members correspond to alphaS variation.
    ## In this case, the errors above are combined PDF+alphaS uncertainties, obtained by
    ## adding in quadrature the PDF and alphaS uncertainties, shown separately below.
    print "xg(PDF) = %.4g + %.4g - %.4g (+- %.4g)" % (unc.central, unc.errplus_pdf, unc.errminus_pdf, unc.errsymm_pdf)
    print "xg(alphaS) = %.4g +- %.4g" % (unc.central, unc.err_as)
print

## Calculate 1-sigma PDF uncertainty on up-quark distribution.
unc = pset.uncertainty(xuAll)
print "xu = %.4g + %.4g - %.4g (+- %.4g)" % (unc.central, unc.errplus, unc.errminus, unc.errsymm)
print "Scale = %.4g" % unc.scale
if pset.errorType.endswith("+as"):
    print "xu(PDF) = %.4g + %.4g - %.4g (+- %.4g)" % (unc.central, unc.errplus_pdf, unc.errminus_pdf, unc.errsymm_pdf)
    print "xu(alphaS) = %.4g +- %.4g" % (unc.central, unc.err_as)
print

## Calculate PDF correlation between gluon and up-quark.
corr = pset.correlation(xgAll, xuAll)
print "Correlation = %.4g" % corr
print

## Generate random values from Hessian best-fit and eigenvector values.
if pset.errorType.startswith("hessian") or pset.errorType.startswith("symmhessian"):
    ## If the ErrorType ends in "+as" exclude the last two members (alphaS variations).
    npdfmem = nmem-2 if pset.errorType.endswith("+as") else nmem
    neigen = npdfmem/2 if pset.errorType.startswith("hessian") else npdfmem
    ## Random numbers are just set to zero here for testing purposes.
    randoms = [0.0 for ir in xrange(neigen)]
    xgrand = pset.randomValueFromHessian(xgAll, randoms)
    xurand = pset.randomValueFromHessian(xuAll, randoms)
    print "Random: xg = %.4g, xu = %.4g" % (xgrand, xurand)

print
