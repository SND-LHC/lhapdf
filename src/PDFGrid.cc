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


  double PDFGrid::_xfxQ2(PID_t id, double x, double q2 ) const {
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


  PDFGrid PDFGrid::load(const string& filepath) {
    ifstream file(filepath.c_str());
    return PDFGrid::load(file);
  }


  PDFGrid PDFGrid::load(istream& file) {
    /// @todo First create a PDFInfo object to parse the top part (and to cascade up to the set info and the global config)

    /// @todo The knots are now part of each subgrid rather than the YAML metadata
    // // Parse grid knots
    // /// @todo Replace with a special vector<double> specialisation of PDF::metadata()?
    // for (YAML::Iterator xsit = head["Xs"].begin(); xsit != head["Xs"].end(); ++xsit) {
    //   double x;
    //   (*xsit) >> x;
    //   grid->_xknots.push_back(x);
    // }
    // /// @todo Replace with a special vector<double> specialisation of PDF::metadata()?
    // for (YAML::Iterator q2sit = head["Q2s"].begin(); q2sit != head["Q2s"].end(); ++q2sit) {
    //   double q2;
    //   (*q2sit) >> q2;
    //   grid->_q2knots.push_back(q2);
    // }

    PDFGrid rtn; //< The return value to-be

    string line;
    int iblock(0), iline(0);
    vector<double> xs, q2s;
    const size_t npid = 11; /// @todo For 5q + g PDFs. Convert to info().flavors().size() once info system exists
    vector< vector<double> > ipid_xfs(npid);
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
        KnotArrayNF& arraynf = rtn._knotarrays[q2s.front()]; //< Reference to newly created subgrid on the return object
        for (size_t ipid = 0; ipid < npid; ++ipid) {
          PID_t pid = ipid; //< @todo Replace with info().flavors()[ipid]; when info() works
          arraynf[pid] = KnotArray1F(xs, q2s); // create the 2D array with the x and Q2 knot positions
          arraynf[pid].xfs().assign(ipid_xfs[ipid].begin(), ipid_xfs[ipid].end()); // populate the xf array
        }
        xs.clear(); q2s.clear(); ipid_xfs.clear();
      }
    }

    // Set default inter/extrapolators
    /// @todo Re-enable when the info system is sorted out
    // const string ipolname = metadata("Interpolator");
    // rtn.setInterpolator(ipolname);
    // const string xpolname = metadata("Extrapolator");
    // rtn.setExtrapolator(xpolname);

    return rtn;
  }


}
