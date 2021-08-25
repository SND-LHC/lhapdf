#! /usr/bin/env python

from __future__ import print_function
import sys

PDF = "CT14nnlo"
PID = 0
if len(sys.argv) > 1:
    PDF = sys.argv[1]
if len(sys.argv) > 2:
    PID = int(sys.argv[2])

import lhapdf
pdfs = lhapdf.mkPDFs(PDF)

import numpy as np
xs = np.append(np.exp(np.linspace(np.log(1e-5), np.log(0.1))), np.linspace(0.1, 1.0))
# print xs

p0 = np.array([pdfs[0].xfxQ(PID, x, 100) for x in xs])
p1 = np.array([pdfs[1].xfxQ(PID, x, 100) for x in xs])
r10 = p1 / p0
print(r10)
# p2 = array([pdfs[2].xfxQ(PID,x,100) for x in xs])
# r20 = p2/p0
# plot(xs, r10, label="EV1/Nom")
# plot(xs, r20, label="EV2/Nom")

from matplotlib import pyplot as plt
plt.figure()
for i in range(1, len(pdfs)): #< skip first intentionally, since computing ratios wrt [0]
    plt.plot(xs, np.array([pdfs[i].xfxQ(PID, x, 100) for x in xs]) / p0, label="EV{:d}/Nom".format(i))
plt.xscale("log")
plt.gca().annotate("PID = {0:d}, Q = 100 GeV".format(PID), (0.2, 0.2), ha="left", xycoords="axes fraction")
# plt.legend()
plt.xlabel("$x$")
plt.ylabel(r"$f(x,Q^2)_\mathrm{{EV}}/f(x,Q^2)_\mathrm{{nom}}$")
# plt.show()
plt.savefig("pdfratios-" + PDF + "-pid" + str(PID) + ".pdf")
