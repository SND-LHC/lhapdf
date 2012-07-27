import numpy
import math

xMin = 0.0
xMax = 1

q2Min = 0.0
q2Max = 1

xN = 20
qN = 20

xs = numpy.linspace( xMin, xMax, xN )
qs = numpy.linspace( q2Min, q2Max, qN )
pids = [0]

print "---"
print " ".join( str(x) for x in xs )
print " ".join( str(q) for q in qs )
print " ".join( str(pid) for pid in pids )

for q in qs:
	for x in xs:
		#print str( 0.1*math.sin(10*x)*math.sin(10*q) )
		print( str( 1.0 ) )
