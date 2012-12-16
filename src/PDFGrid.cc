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


  void PDFGrid::setInterpolator(Interpolator* ipol) {
    _interpolator = ipol;
    _interpolator->bind(*this);
  }

  void PDFGrid::setDefaultInterpolator() {
    const string ipolname = metadata("Interpolator");
    setInterpolator(createInterpolator(ipolname));
  }

  void PDFGrid::setExtrapolator(Extrapolator* xpol) {
    _extrapolator = xpol;
    _extrapolator->bind(*this);
  }

  void PDFGrid::setDefaultExtrapolator() {
    const string xpolname = metadata("Extrapolator");
    setExtrapolator(createExtrapolator(xpolname));
  }


  double PDFGrid::xfxQ2(PID_t id, double x, double q2 ) const {
    /// @todo Move these first three checks into PDF: should *always* be called
    // Physical x range check
    if (!inPhysicalRangeX(x)) {
      stringstream error;
      error << "Unphysical x range given: " << x;
      throw runtime_error( error.str() );
    }
    // Physical Q2 range check
    if (!inPhysicalRangeQ2(q2)) {
      stringstream error;
      error << "Unphysical Q2 range given: " << q2;
      throw runtime_error( error.str() );
    }
    // Undefined PIDs
    if (!hasFlavor(id)) {
      stringstream error;
      error << "Undefined flavour requested: " << id;
      throw runtime_error( error.str() );
    }

    // Decide whether to use interpolation or extrapolation
    if (inRangeXQ2(x, q2)) {
      if (!hasInterpolator()) throw runtime_error("Undefined interpolator");
      return _interpolator->interpolateQ2(*this, id, x, q2);
    } else {
      if (!hasExtrapolator()) throw runtime_error("Undefined extrapolator");
      return _extrapolator->extrapolateQ2(*this, id, x, q2);
    }
  }



  inline void tokenize(const string& line, AxisKnots& knots) {
    size_t start = 0;
    size_t pos = 0;
    do {
      pos = line.find_first_of( " ", start );
      knots.push_back( atof( line.substr( start, (pos - start) ).c_str() ) );
      start = pos+1;
    } while (pos != string::npos);
  }


  PDFGrid* PDFGrid::load( PDFGrid* grid, const YAML::Node& head, ifstream& file ) {
    // Parse grid knots
    for (YAML::Iterator xsit = head["Xs"].begin(); xsit != head["Xs"].end(); ++xsit) {
      double x;
      (*xsit) >> x;
      grid->_xknots.push_back(x);
    }

    for (YAML::Iterator q2sit = head["Q2s"].begin(); q2sit != head["Q2s"].end(); ++q2sit) {
      double q2;
      (*q2sit) >> q2;
      grid->_q2knots.push_back(q2);
    }

    // Parse grid data
    // Allocate flavor data
    vector<PID_t>::const_iterator piditer = grid->_set->flavors().begin();
    for( ; piditer != grid->_set->flavors().end(); ++piditer ) {
      grid->flavors[*piditer] = new double[grid->xKnots().size() * grid->q2Knots().size()];
    }

    // Parse grid lines
    unsigned int cline = 0;
    for(; cline < grid->xKnots().size()*grid->q2Knots().size(); ++cline ) {
      // Read a line
      if( !file.good() ) {
        stringstream error;
        error << "ifstream ran out of data @ " << cline;
        throw runtime_error( error.str() );
      }
      string line;
      getline( file, line );

      // Parsing individual grid line
      const char* cstr = line.c_str();
      char str[line.size()+1];  /// @todo Not ISO C++
      memcpy( str, cstr, (line.size()+1)*sizeof(char) );

      char *token, *prog;
      unsigned int flavor=0;
      token = strtok_r( str, " ", &prog );
      while( token != NULL ) {
        // Process token
        grid->flavors[grid->_set->flavors()[flavor]][cline] = atof( token );
        token = strtok_r( NULL, " ", &prog );
        ++flavor;
      }
    }

    // Set default inter/extra/polators
    grid->setDefaultInterpolator();
    grid->setDefaultExtrapolator();

    return grid;
  }
}
