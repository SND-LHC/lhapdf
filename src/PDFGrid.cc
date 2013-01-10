#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/Factories.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstring>

using namespace std;

namespace LHAPDF {


  double PDFGrid::_xfxQ2(int id, double x, double q2 ) const {
    // Decide whether to use interpolation or extrapolation... the sanity checks
    // are done in the public PDF::xfxQ2 function.
    if (inRangeXQ2(x, q2)) {
      if (interpolator() == 0) throw GridError("Undefined interpolator");
      return interpolator()->interpolateXQ2(id, x, q2);
    } else {
      if (extrapolator() == 0) throw GridError("Undefined extrapolator");
      return extrapolator()->extrapolateXQ2(id, x, q2);
    }
  }


  void PDFGrid::_loadData(const string& mempath) {
    string line;
    int iblock(0), iline(0);
    vector<double> xs, q2s;
    const size_t npid = info().metadata< vector<int> >("Flavors").size(); //< @todo Convert to PDF::flavors().size() once it exists
    vector< vector<double> > ipid_xfs(npid);

    try {
      ifstream file(mempath.c_str());
      while (getline(file, line)) {
        iline += 1;
        if (iblock > 0) { // Block 0 is the metadata, which we ignore here
          double token;
          istringstream tokens(line);
          if (iline == 1) { // x knots line
            while (tokens >> token) xs.push_back(token);
          } if (iline == 2) { // Q2 knots line
            while (tokens >> token) q2s.push_back(token);
          } else {
            if (iline == 3) { // on the first line of the xf block, resize the arrays
              for (size_t ipid = 0; ipid < npid; ++ipid) { ipid_xfs[ipid].reserve(xs.size() * q2s.size()); }
            }
            int ipid = 0;
            while (tokens >> token) {
              ipid_xfs[ipid].push_back(token);
              ipid += 1;
            }
          }
        } else if (line == "---") { // This is the block divider line
          iblock += 1;
          iline = 0;
          KnotArrayNF& arraynf = _knotarrays[q2s.front()]; //< Reference to newly created subgrid on the return object
          for (size_t ipid = 0; ipid < npid; ++ipid) {
            int pid = ipid; //< @todo Replace with info().flavors()[ipid]; when info() works
            arraynf[pid] = KnotArray1F(xs, q2s); // create the 2D array with the x and Q2 knot positions
            arraynf[pid].xfs().assign(ipid_xfs[ipid].begin(), ipid_xfs[ipid].end()); // populate the xf array
          }
          xs.clear(); q2s.clear(); ipid_xfs.clear();
        }
      }
    } catch (std::exception& e) {
      throw ReadError("Read error while parsing " + mempath + " as a PDFGrid data file");
    }

    // Set default inter/extrapolators
    const string ipolname = info().metadata("Interpolator");
    setInterpolator(ipolname);
    const string xpolname = info().metadata("Extrapolator");
    setExtrapolator(xpolname);
  }


}
