#! /usr/bin/env python

import sys
PDF = "CT14nnlo"
PID = 0
if len(sys.argv) > 1:
    PDF = sys.argv[1]
if len(sys.argv) > 2:
    PID = int(sys.argv[2])

import lhapdf
pdfs = lhapdf.mkPDFs(PDF)

import pylab as pl
xs = pl.append(pl.exp(pl.linspace(pl.log(1e-5), pl.log(0.1))), pl.linspace(0.1,1.0))
#print xs

p0 = pl.array([pdfs[0].xfxQ(PID,x,100) for x in xs])
p1 = pl.array([pdfs[1].xfxQ(PID,x,100) for x in xs])
r10 = p1/p0
print r10
# p2 = array([pdfs[2].xfxQ(PID,x,100) for x in xs])
# r20 = p2/p0
# plot(xs, r10, label="EV1/Nom")
# plot(xs, r20, label="EV2/Nom")

for i in range(1, len(pdfs)):
    pl.plot(xs, pl.array([pdfs[i].xfxQ(PID,x,100) for x in xs])/p0, label="EV{:d}/Nom".format(i))
pl.xscale("log")
pl.gca().annotate("PID = {0:d}, Q = 100 GeV".format(PID), (0.2,0.2), ha="left", xycoords="axes fraction")
#pl.legend()
pl.xlabel("$x$")
pl.ylabel(r"$f(x,Q^2)_\mathrm{{EV}}/f(x,Q^2)_\mathrm{{nom}}$")
#pl.show()
pl.savefig("pdfratios-" + PDF + "-pid" + str(PID) + ".pdf")
