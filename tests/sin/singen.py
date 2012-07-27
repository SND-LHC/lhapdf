import numpy
import math

xN = 16
qN = 16

xMin = 0.0
xMax = 1

q2Min = 0.0
q2Max = 1

xs = numpy.linspace( xMin, xMax, xN )
qs = numpy.linspace( q2Min, q2Max, qN )
pids = [0]

print "---"
print " ".join( str(x) for x in xs )
print " ".join( str(q) for q in qs )
print " ".join( str(pid) for pid in pids )

for q in qs:
	for x in xs:
		print str( 0.04*math.sin(40*x)*math.sin(40*q) )
		#print str( 0.1*math.sin(10*x) )
		#print str( 0.1*math.sin(10*q) )
		#print( str( 1.0 ) )
