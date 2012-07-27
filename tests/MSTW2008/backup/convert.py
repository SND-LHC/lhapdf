import lhapdf
import numpy

#Converting MSTW2008lo90cl set
lhapdf.initPDFSetByName("MSTW2008lo90cl.LHgrid")
#Use only central PDF
lhapdf.initPDF(0)

#Define X Points
#xpoints = [ 1e-6,	2e-6,	4e-6,	6e-6,	8e-6,
#			1e-5,	2e-5,	4e-5,	6e-5,	8e-5,
#			1e-4,	2e-4,	4e-4,	6e-4,	8e-4,
#			1e-3,	2e-3,	4e-3,	6e-3,	8e-3,
#			1e-2,	1.4e-2,	2e-2,	3e-2,	4e-2,
#			6e-2,	8e-2,	0.1,	0.125,	0.15,
#			0.175,	0.2,	0.225,	0.25,	0.275,
#			0.3,	0.325,	0.35,	0.375,	0.4,
#			0.425,	0.45,	0.475,	0.5,	0.525,
#			0.55,	0.575,	0.6,	0.625,	0.65,
#			0.675,	0.7,	0.725,	0.75,	0.775,
#			0.8,	0.825,	0.85,	0.875,	0.9,
#			0.925,	0.95,	0.975,	1.0 ]

#Define Q Points
#qpoints = [ 1.0,	1.25,	1.5,	0.0,	0.0,
#			2.5,	3.2,	4,		5,		6.4,
#			8.0, 	1.0e1,	1.2e1,	0.0,	0.0,
#			2.6e1,	4.0e1,	6.4e1,	1.0e2,	1.6e2,
#			2.4e2,	4.0e2,	6.4e2, 	1.0e3,	1.8e3,
#			3.2e3,	5.6e3,	1.0e4,	1.8e4,	3.2e4,
#			5.6e4, 1.0e5,	1.8e5,	3.2e5,	5.6e5,
#			1.0e6,	1.8e6,	3.2e6,	5.6e6, 	1.0e7,
#			1.8e7,	3.2e7,	5.6e7,	1.0e8,	1.8e8,
#			3.2e8,	5.6e8, 1.0e9 ]

#Define PID translation
pids = [ -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 ]

#Generate our own grid locations, because the others are not well defined yet
xpoints = numpy.linspace(1e-6, 1.0, 64)
q2points = numpy.linspace(1.0, 1.0e9, 48)

xfxs = []
for x in xpoints:
        for q2 in q2points:
                xfxs.append([lhapdf.xfx(x, numpy.sqrt(q2), k) for k in xrange(-6,  7)])

f = open('0.LHm', 'w')

f.write(" ".join(str(x) for x in xpoints))
f.write('\n')

f.write(" ".join(str(q2) for q2 in q2points))
f.write('\n')

f.write(" ".join(str(pid) for pid in pids))
f.write('\n')

for line in xfxs:
	print >>f, " ".join(str(i) for i in line)
			
f.close()