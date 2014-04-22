#! /usr/bin/env python

import lhapdf
lhapdf.pathsAppend(".")

#p_lh = lhapdf.mkPDF("ATLAS-epWZ12-EIG", 0)
p_lh = lhapdf.mkPDF("ATLAS-epWZ12-EIG_LH", 0)
p_hf = lhapdf.mkPDF("ATLAS-epWZ12-EIG_HF", 0)

import numpy
q = 1.9
xs = list(numpy.logspace(-3, -0.01, 100))
xs = list(numpy.linspace(0.01, 0.99, 98))

xfgs_lh = [p_lh.xfxQ2(21, x, q) for x in xs]
xfgs_hf = [p_hf.xfxQ2(21, x, q) for x in xs]

def fmt(xs):
    return ["%4g" % x for x in xs]

print
print fmt(xs), "\n"
print fmt(xfgs_lh), "\n"
print fmt(xfgs_hf), "\n"
print fmt([xfgs_hf[i]/xfgs_lh[i] - 1.0 for i in xrange(len(xs))]), "\n"
