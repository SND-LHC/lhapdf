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
    /// Decide whether to use interpolation or extrapolation... the sanity checks
    /// are done in the public PDF::xfxQ2 function.
     //cout << "From GridPDF[0]: x = " << x << ", Q2 = " << q2 << endl;
    if (inRangeXQ2(x, q2)) {
       //cout << "From GridPDF[ipol]: x = " << ", Q2 = " << q2 << x << endl;
      return interpolator().interpolateXQ2(id, x, q2);
    } else {
       //cout << "From GridPDF[xpol]: x = " << x << ", Q2 = " << q2 << endl;
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
        // Trim the current line to ensure that there is no effect of leading spaces, etc.
        trim(line);

        // If the line is commented out, increment the line number but not the block line
        iline += 1;
        if (line.find("#") == 0) continue;
        iblockline += 1;

        if (line != "---") { // if we are not on a block separator line...

          // Block 0 is the metadata, which we ignore here
          if (iblock == 0) continue;

          // Debug
          //cout << line << " @ " << iline << " = block line #" << iblockline << endl;

          // Parse the data lines
          double token;
          istringstream tokens(line);
          if (iblockline == 1) { // x knots line
            while (tokens >> token) xs.push_back(token);
          } else if (iblockline == 2) { // Q2 knots line
            while (tokens >> token) q2s.push_back(token);
          } else if (iblockline == 3) { // internal flavor IDs line
            // DO NOTHING FOR NOW: only handling this for prospective forward compatibility
            /// @todo Handle internal partial flavour representations
            //while (tokens >> token) intflavors.push_back(token);
          } else {
            if (iblockline == 4) { // on the first line of the xf block, resize the arrays
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
          if (iblock > 0 && iblockline - 1 != int(xs.size()*q2s.size()) + 3)
            throw ReadError("PDF grid data error on line " + to_str(iline) + ": " +
                            to_str(iblockline-1) + " data lines were seen in block " + to_str(iblock-1) +
                            " but " + to_str(xs.size()*q2s.size() + 3) + " expected");

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

          // Register data from the previous (>0th) block into the GridPDF data structure
          // KnotArrayNF arraynf;
          KnotArrayNF& arraynf = _knotarrays[q2s.front()]; //< Reference to newly created subgrid object
          for (size_t ipid = 0; ipid < npid; ++ipid) {
            const int pid = flavors()[ipid];
            arraynf[pid] = KnotArray1F(xs, q2s); // create the 2D array with the x and Q2 knot positions
            // Populate the xf data array (version 1: neat but portability not guaranteed?)
            arraynf[pid].xfs().assign(ipid_xfs[ipid].begin(), ipid_xfs[ipid].end());
            // // Populate the xf data array (version 2: one loop across list, converting to 2 subindices)
            // for (size_t i = 0; i < ipid_xfs[ipid].size(); ++i) {
            //   arraynf[pid].xfs()[i / q2s.size()][i % q2s.size()] = ipid_xfs[ipid][i];
            // }
            // // Populate the xf data array (version 3: two loops across subindices, constructing list index)
            // for (ix = 0; ix < xs.size(); ++ix) {
            //   for (iq = 0; iq < q2s.size(); ++iq) {
            //     arraynf[pid].xfs()[ix][iq] = ipid_xfs[ipid][iq*q2s.size() + ix];
            //   }
            // }
          }
          //cout << _knotarrays.size() << endl;
          xs.clear(); q2s.clear();
          for (size_t ipid = 0; ipid < npid; ++ipid) ipid_xfs[ipid].clear();

        }
      }
    } catch (Exception& e) {
      throw;
    } catch (std::exception& e) {
      throw ReadError("Read error while parsing " + mempath.native() + " as a GridPDF data file");
    }

  }


}
