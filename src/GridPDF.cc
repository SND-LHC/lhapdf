#include "LHAPDF/GridPDF.h"
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


  double GridPDF::_xfxQ2(int id, double x, double q2 ) const {
    // Decide whether to use interpolation or extrapolation... the sanity checks
    // are done in the public PDF::xfxQ2 function.
    if (inRangeXQ2(x, q2)) {
      return interpolator().interpolateXQ2(id, x, q2);
    } else {
      return extrapolator().extrapolateXQ2(id, x, q2);
    }
  }


  void GridPDF::_loadData(const path& mempath) {
    string line;
    int iblock(0), iblockline(0), iline(0);
    vector<double> xs, q2s;
    const size_t npid = flavors().size();
    vector< vector<double> > ipid_xfs(npid);

    try {
      ifstream file(mempath.c_str());
      while (getline(file, line)) {
        iline += 1;
        iblockline += 1;

        // Trim the current line to ensure that there is no effect of leading spaces, etc.
        trim(line);

        if (line != "---") { // if we are not on a block separator line...

          // Block 0 is the metadata, which we ignore here
          if (iblock == 0) continue;

          // Parse the data lines
          double token;
          istringstream tokens(line);
          if (iblockline == 1) { // x knots line
            while (tokens >> token) xs.push_back(token);
          } else if (iblockline == 2) { // Q2 knots line
            while (tokens >> token) q2s.push_back(token);
          } else {
            if (iblockline == 3) { // on the first line of the xf block, resize the arrays
              for (size_t ipid = 0; ipid < npid; ++ipid) { ipid_xfs[ipid].reserve(xs.size() * q2s.size()); }
            }
            size_t ipid = 0;
            while (tokens >> token) {
              ipid_xfs[ipid].push_back(token);
              ipid += 1;
            }
            // Check that each line has many tokens as there should be flavours
            if (ipid != npid)
              throw ReadError("PDF grid data error on line " + to_str(iline) + ": " + to_str(ipid) +
                              " flavor entries seen but " + to_str(npid) + " expected");
          }

        } else { // we *are* on a block separator line

          // Check that the expected number of data lines were seen in the last block
          if (iblock > 0 && iblockline - 1 != int(xs.size()*q2s.size()) + 2)
            throw ReadError("PDF grid data error on line " + to_str(iline) + ": " +
                            to_str(iblockline-1) + " data lines were seen in block " + to_str(iblock-1) +
                            " but " + to_str(xs.size()*q2s.size() + 2) + " expected");

          // Increment/reset the block and line counters
          iblock += 1;
          iblockline = 0;

          // Escape here if we've just finished reading the 0th (metadata) block
          if (iblock == 1) continue;

          // Die with an assert if the block was of zero size
          /// @todo Convert to throwing some exception? Is this ever allowable?
          assert(xs.size() > 0);
          assert(q2s.size() > 0);
          assert(ipid_xfs.size() > 0);

          /// @todo Define the ordering of the values in x and Q2 indexing

          // Register data from the previous (>0th) block into the GridPDF data structure
          // KnotArrayNF arraynf;
          KnotArrayNF& arraynf = _knotarrays[q2s.front()]; //< Reference to newly created subgrid object
          for (size_t ipid = 0; ipid < npid; ++ipid) {
            int pid = flavors()[ipid];
            arraynf[pid] = KnotArray1F(xs, q2s); // create the 2D array with the x and Q2 knot positions
            arraynf[pid].xfs().assign(ipid_xfs[ipid].begin(), ipid_xfs[ipid].end()); // populate the xf array
          }
          cout << q2s.size() << endl;
          cout << q2s.front() << endl;
          // _knotarrays[q2s.front()] = arraynf; //< @todo Prefer getting the reference as above, to avoid re-copying a big data array
          xs.clear(); q2s.clear(); ipid_xfs.clear();

        }
      }
    } catch (Exception& e) {
      throw;
    } catch (std::exception& e) {
      throw ReadError("Read error while parsing " + mempath.native() + " as a GridPDF data file");
    }

  }


}
