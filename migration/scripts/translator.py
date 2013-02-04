import lhapdf
import numpy

lhapdf.initPDFSetByName("MSTW2008nnlo90cl.LHgrid")
lhapdf.initPDF(0)

nMembers = lhapdf.numberPDF()

xMin = lhapdf.getXmin(0)
xMax = lhapdf.getXmax(0)

q2Min = lhapdf.getQ2min(0)
q2Max = lhapdf.getQ2max(0)

"""
for i in nMembers:
	xMinTemp = lhapdf.getXmin(i)
	xMaxTemp = lhapdf.getXmax(i)
	q2MinTemp = lhapdf.getQ2min(i)
	q2MaxTemp = lhapdf.getQ2max(i)
	if xMinTemp < xMin:
		xMin = xMinTemp
	if xMaxTemp > xMax:
		xMax = xMaxTemp
	if q2MinTemp < q2Min:
		q2Min = q2MinTemp
	if q2MaxTemp < q2Max:
		q2Max = q2MaxTemp
"""

nx = 100
nq2 = 100

xStep = (xMax - xMin)/nx
q2Step = (q2Max - q2Min)/nq2

print xStep
print q2Step

xfxs = []

x = xMin
q2 = q2Min
pids = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]	

xs = numpy.logspace(numpy.log(xMin), numpy.log(xMax), nx)
qs = numpy.linspace(q2Min, q2Max, nq2)

for x in xs:
	for q in qs:
		xfxs.append([lhapdf.xfx(x, numpy.sqrt(q), k) for k in xrange(-6,  7)])

"""
for i in xfxs:
	print i
"""

f = open('member_0.LHm', 'w')

f.write(" ".join(str(x) for x in xs))
f.write('\n')

f.write(" ".join(str(q) for q in qs))
f.write('\n')

f.write(" ".join(str(pid) for pid in pids))
f.write('\n')

for line in xfxs:
	print >>f, " ".join(str(i) for i in line)