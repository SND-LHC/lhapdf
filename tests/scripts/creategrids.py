import gridhack
import lhapdf
import numpy
import os

def setup():
    lhapdf.initPDFSetByName('cteq66.LHgrid')
    return gridhack.get_grid()

def xfx_values(lhapids, xs, q2s):
    xfxs = []
    for q2 in q2s:
        for x in xs:
            xfxs.append([lhapdf.xfx(x, q2, pid) for pid in lhapids])
    return xfxs

def active_flavours(xfxs):
    pids = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
    xfxs_array = numpy.array(xfxs)
    pids = [pid for i, pid in enumerate(pids) if xfxs_array[:,i].any()]
    return pids


# Select PDF sets to translates here
pdf_sets = ['MSTW2008lo90cl_nf3', 'MSTW2008lo90cl', 'NNPDF23_nlo_as_0119', 'CT10', 'cteq66', 'cteq61']

# Declare the base path directory to LHGrid files
dir = '/phys/linux/s0821167/Summer Project/Interpolator/lhapdfv6/tests/'

# Declare meta filename, the same for every set.
meta_filename = 'meta.LHinfo'

# For every PDF set
for pdf_set in pdf_sets:

    # Initialise PIDs 
    lhapids = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]

    # Initialise LHgrid file
    input = pdf_set + '.LHgrid'
    lhapdf.initPDFSetByName(input)

    # Determine correct output directory for set
    output_dir = dir + pdf_set + '/'
    
    # Create directory
    try:
        os.mkdir(output_dir)
    except OSError:
        pass # already exists

    # Use gridhack python wrapper to call LHAPDF for exact x and q2 grid points
    xs, q2s = gridhack.get_grid()

    # Determine active flavours
    xfxs = xfx_values(lhapids, xs, q2s)
    lhapids = active_flavours(xfxs)        

    # Create meta file for set
    with open(os.path.join(output_dir, meta_filename), 'wb') as meta_output:
        print >>meta_output, '---'
        print >>meta_output, 'SetName: ' + pdf_set
        print >>meta_output, 'SetDescription: '
        print >>meta_output
        print >>meta_output, 'SetLambdaQCD: '
        print >>meta_output, 'SetOrder: '
        print >>meta_output
        print >>meta_output, 'SetFlavours: '
        print >>meta_output
        print >>meta_output, 'SetInterpolator: '
        print >>meta_output, 'SetExtrapolator: '
        print >>meta_output, '---'

    # Find the number of pdf members in pdf set
    n_members = lhapdf.numberPDF()

    for member in range(n_members):

        # Load member
        lhapdf.initPDF(member)

        # Set up new output grid format file
        output_file = 'mbr_' + str(member)+ '.LHgm'

        # Set member name to be central if it is the 0th member, errors
        # otherwise
        member_name = ''
        if member == 0:
            member_name = 'Central'
        else:
            member_name = 'Errors'

        # Get xfx values
        xfxs = xfx_values(lhapids, xs, q2s)

        # Print results to file:
        # ---------------------
        # xs
        # q2s
        # xfxs
        # ---------------------
        with open(os.path.join(output_dir, output_file), 'wb') as output:
            print >>output, 'MemberName: ' + member_name
            print >>output, 'MemberID: ' + str(member)
            print >>output, 'Xs: ' + str(xs)
            print >>output, 'Q2s: ' + str(q2s)           
            print >>output, '---'
            for line in xfxs:
                print >>output, " ".join(str(i) for i in line)
