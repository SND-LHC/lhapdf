#! /usr/bin/env python

import lhapdf, numpy
import os, sys, optparse

def getxfs(lhapids, xs, qs):
    xfs = []
    for x in xs:
        for q in qs:
            xfs.append([lhapdf.xfx(x, q, pid) for pid in lhapids])
    return xfs

def active_flavours(xfxs):
    pids = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
    xfxs_array = numpy.array(xfxs)
    pids = [pid for i, pid in enumerate(pids) if xfxs_array[:,i].any()]
    return pids


## Parse the command line to get the set names for migration (in LHAPDF5 format)
# e.g. 'MSTW2008lo90cl_nf3.LHgrid', 'MSTW2008lo90cl.LHgrid', 'NNPDF23_nlo_as_0119.LHgrid', 'CT10.LHgrid', 'cteq66.LHgrid', 'cteq6ll.LHpdf', 'cteq6ll.LHgrid', ...
p = optparse.OptionParser(usage="%prog <setname> [<setname2> ...]")
opts, args = p.parse_args()


# For every PDF set
for lha5name in args:
    setname = os.path.splitext(lha5name)[0]
    print "Migrating %s -> %s" % (lha5name, setname)

    ## Create output dir if needed
    if not os.path.exists(setname):
        os.mkdir(setname)

    ## Initialize LHAPDF for this set
    lhapdf.initPDFSetByName(lha5name)
    # TODO: Get alpha_s anchors here (for grid Qs?) to write into set metadata

    ## Create meta file for set
    metapath = os.path.join(setname, setname + '.info')
    f = open(metapath, "w")
    f.write("SetDesc: \n")
    f.write("NumMembers: %d\n" % lhapdf.numberPDF())
    f.write("Flavors: [-5,-4,-3,-2,-1,1,2,3,4,5,21]\n")
    f.write("ErrorType: \n")
    f.write("Format: lhagrid1\n")
    f.write("AlphaS_MZ: \n")
    f.write("LambdaQCD: \n")
    f.write("OrderQcd: \n")
    f.close()

    ## Iterate over each member in the set
    for member in xrange(lhapdf.numberPDF()):
        lhapdf.initPDF(member)

        ## Get the x and Q anchor point arrays
        try:
            ## Use LHAPDF5 to get the exact x and q2 grid points (requires LHAPDF 5.9)
            import gridhack
            xs, q2s = gridhack.get_grid()
        except:
            xs = numpy.logspace(-8, -0.001, 30)
            q2s = numpy.logspace(0.1, 8, 50)
        qs = numpy.sqrt(numpy.array(q2s))

        ## Determine the active (non-zero) flavours in the grid (will be -5..5 or -6..6)
        ## NB. 0 is listed last since it will be translated to 21 in the PDG scheme
        #lhapids = active_flavours(xfxs)
        #lhapids = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 0]
        lhapids = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 0]

        ## Get the xf values for this PDF member
        xfs = getxfs(lhapids, xs, qs)

        ## Work out the member file path and open it for writing
        memname = setname + ('_%04d.lha' % member)
        mempath = os.path.join(setname, memname)
        f = open(mempath, "w")
        f.write("PdfDesc: \n")
        f.write("PdfType: %s\n" % ("central" if member == 0 else "error"))
        f.write("---\n")
        ## Write x points
        line = " ".join("%2.6e" % x for x in xs)
        f.write(line + "\n")
        ## Write Q2 points
        line = " ".join("%2.6e" % q2 for q2 in q2s)
        f.write(line + "\n")
        ## Write block of xf values
        for xfs_xq in xfs:
            line = ""
            for xf in xfs_xq:
                if xf == 0.0: xf = 0.0 # remove occurences of negative zero
                line += "%2.6e " % xf
            f.write(line.strip() + "\n")

        f.write("---\n")
        f.close()
