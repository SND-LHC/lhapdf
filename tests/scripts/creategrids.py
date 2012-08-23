import gridhack
import lhapdf
import numpy
import os
import os.path

# Select PDF sets to translates here
pdf_sets = ['MSTW2008lo90cl_nf3']

base_dir = "."

# For every PDF set
for pdf in pdf_sets:
	# Initialise LHgrid file
	input = pdf
	input += '.LHgrid'
	lhapdf.initPDFSetByName(input)
	
	setpath = os.path.join(base_dir,pdf)
	try:
		os.mkdir(setpath)
	except:
		pass

	for mem in range(lhapdf.numberPDF()):
		# Set up new output grid format file
		output = os.path.join(setpath, "mbr_%s.LHm" % mem)

		# Use gridhack python wrapper to call LHAPDF for exact x and q2 grid points
		xs, q2s = gridhack.get_grid()

		# Hard-coded pids, excluding none existent ttbar
		pids = [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]
		lhapids = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]

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
			print >>output, "MemberType: PDFGrid"
			print >>output, "MemberName: err%s" % mem 
			print >>output, "MemberID: %s" % mem
			print >>output, "Xs: %r" % xs
			print >>output, "Q2s: %r" % q2s
			
			print >>output, "---"

			for line in xfxs:
				print >>output, " ".join(str(i) for i in line)
