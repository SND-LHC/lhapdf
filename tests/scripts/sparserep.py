from sys import stdin
import optparse

parser = optparse.OptionParser()
parser.add_option("-d", "--delimiter", dest="delimiter", help="Specify a delimiter for output", default=" ")

opts, args = parser.parse_args()

try:
	filename = args[0]
except IndexError, e:
 	raise RuntimeError('No filename')

lines = open(filename)

xs = map( float, next(lines).split() )
qs = map( float, next(lines).split() )
pids = map( int, next(lines).split() )

from itertools import product, izip

files = []

for (q, x), line in izip( product( qs, xs ), lines ):
	values = map( float, line.split() )
	print x, q, opts.delimiter.join(map(str, values))

for file in files:
	file.close()
lines.close()
