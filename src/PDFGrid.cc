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



  inline void tokenize(const string& line, vector<double>& knots) {
    size_t start = 0;
    size_t pos = 0;
    do {
      pos = line.find_first_of( " ", start );
      knots.push_back( atof( line.substr( start, (pos - start) ).c_str() ) );
      start = pos+1;
    } while (pos != string::npos);
  }



  /// @todo Should just be a filename as argument: why expose YAML or pass a PDFGrid?
  PDFGrid* PDFGrid::load( PDFGrid* grid, const YAML::Node& head, ifstream& file) {
    // Parse grid knots
    /// @todo Replace with a special vector<double> specialisation of PDF::metadata()?
    for (YAML::Iterator xsit = head["Xs"].begin(); xsit != head["Xs"].end(); ++xsit) {
      double x;
      (*xsit) >> x;
      grid->_xknots.push_back(x);
    }
    /// @todo Replace with a special vector<double> specialisation of PDF::metadata()?
    for (YAML::Iterator q2sit = head["Q2s"].begin(); q2sit != head["Q2s"].end(); ++q2sit) {
      double q2;
      (*q2sit) >> q2;
      grid->_q2knots.push_back(q2);
    }

    // Parse grid data
    // Allocate knot arrays

    for (vector<PID_t>::const_iterator fl = grid->flavors().begin(); fl != grid->flavors().end(); ++fl) {
      /// @todo Remember to call delete[] when destructing...
      grid->_ptdata[*fl] = new double[grid->xKnots().size() * grid->q2Knots().size()];
    }

    // Parse grid lines
    size_t cline = 0;
    for (; cline < grid->xKnots().size()*grid->q2Knots().size(); ++cline) {
      // Read a line
      if (!file.good()) {
        stringstream error;
        error << "ifstream ran out of data @ " << cline;
        throw ReadError(error.str());
      }
      string line;
      getline(file, line);

      // Parsing individual grid line
      const char* cstr = line.c_str();
      char str[line.size()+1];  /// @todo Not ISO C++. Fix!
      memcpy( str, cstr, (line.size()+1)*sizeof(char) );

      char *token, *prog;
      unsigned int flavor = 0;
      token = strtok_r( str, " ", &prog );
      while (token != NULL) {
        // Process token
        /// @todo Massively fix and improve...
        // AB tmp disable: grid->_ptdata[grid->_set->flavors()[flavor]][cline] = atof( token );
        token = strtok_r(NULL, " ", &prog);
        ++flavor;
      }
    }

    // Set default inter/extrapolators
    /// @todo Re-enable when the info system is sorted out
    // const string ipolname = metadata("Interpolator");
    // grid->setInterpolator(ipolname);
    // const string xpolname = metadata("Extrapolator");
    // grid->setExtrapolator(xpolname);

    return grid;
  }
}
