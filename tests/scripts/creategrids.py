import gridhack
import lhapdf
import numpy

# Select PDF sets to translates here
pdf_sets = ['MSTW2008lo90cl', 'NNPDF23_nlo_as_0119', 'CT10', 'cteq66', 'cteq61']

# For every PDF set
for pdf in pdf_sets:

	# Initialise LHgrid file
	input = pdf
	input += '.LHgrid'
	lhapdf.initPDFSetByName(input)

	# Set up new output grid format file
	output = pdf
	output += '.LHm'

	# Use gridhack python wrapper to call LHAPDF for exact x and q2 grid points
	xs, q2s = gridhack.get_grid()

	# Hard-coded pids, excluding none existent ttbar
	pids = [-5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5]
	lhapids = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]

	# Use C++ LHAPDF wrapper to return xfx values at our exact grid points
	xfxs = []
	for q2 in q2s:
		for x in xs:
			xfxs.append([lhapdf.xfx(x, sqrt(q2), pid) for pid in lhapids])

	# Print results to file:
	# ---------------------
	# xs
	# q2s
	# pids
	# xfxs
	# ---------------------
	with open(output, 'w') as output:
		print >>output, " ".join(str(x) for x in xs)
		print >>output, " ".join(str(q2) for q2 in q2s)
		print >>output, " ".join(str(pid) for pid in pids)
		for line in xfxs:
			print >>output, " ".join(str(i) for i in line)
