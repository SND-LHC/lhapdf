import lhapdf
import numpy

lhapdf.initPDFSetByName("MRST2001lo.LHgrid")
lhapdf.initPDF(0)

xMin = 1e-10
xMax = lhapdf.getXmax(0)

qMin = 1e-10
#qMin = 0
#qMax = numpy.sqrt(lhapdf.getQ2max(0))
qMax = 1e10

nx = 100	
nq = 100

pids = [-5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5]	

xlogs = numpy.linspace(numpy.log10(xMin), numpy.log10(xMax), nx)
qlogs = numpy.linspace(numpy.log10(qMin), numpy.log10(qMax), nq)

xlins = 10**xlogs
qlins = 10**qlogs

xfxs = []
for q in qlins:
	for x in xlins:
		xfxs.append([lhapdf.xfx(x, q, pid) for pid in xrange(-5, 6)])

print " ".join(str(x) for x in xlogs)

print " ".join(str(q) for q in qlogs)

print " ".join(str(pid) for pid in pids)

for line in xfxs:
	print " ".join(str(i) for i in line)
