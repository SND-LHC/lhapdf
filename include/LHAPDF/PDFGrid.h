#pragma once

#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Types.h"
#include "LHAPDF/Factories.h"
#include <vector>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <exception>
#include <fstream>
#include "yaml-cpp/yaml.h"

namespace LHAPDF {


  // Forward declarations
  class Interpolator;
  class Extrapolator;


  /// @brief A PDF defined via an interpolation grid
  class PDFGrid : public PDF {
  public:


    /// Constructor
    PDFGrid(PDFSet* setp = 0)
      : //_set(setp), /// @todo Why didn't this work?!
        _interpolator(0) _extrapolator(0)
    {
      /// @todo Parse metadata file, create set if needed, set up alpha_s object, etc.
    }


    /// Destructor
    ~PDFGrid() {
      delete _interpolator;
      delete _extrapolator;
    }


    /// @brief PDF xf(x,Q2) value
    double xfxQ2( PID_t, double x, double q2 ) const;


    /// Metadata
    //@{

    /// Get the list of available flavours by PDG ID code.
    /// @todo Or get the flavour list from the set?
    std::vector<PID_t> flavors() const {
      std::vector<PID_t> rtn;
      for (std::map<PID_t, double*>::const_iterator i = _rawdata.begin(); i != _rawdata.end(); ++i) {
        rtn.push_back(i->first);
      }
      return rtn;
    }

    /// Check if x is in the grid range
    bool inRangeX(double x) const {
      if (x < xKnots().front()) return false;
      if (x > xKnots().back()) return false;
      return true;
    }

    /// Check if q2 is in the grid range
    bool inRangeQ2(double q2) const {
      if (q2 < q2Knots().front()) return false;
      if (q2 > q2Knots().back()) return false;
      return true;
    }

    //@}


    /// @name Interpolators and extrapolators
    //@{

    /// Sets the Interpolator to be used for interpolating between grid knots.
    void setInterpolator(Interpolator* ipol);

    /// Sets the default interpolator.
    /// @todo Needed on public API? Just read from grid file.
    void setDefaultInterpolator();

    /// Sets the Extrapolator to be used for extrapolating outside the grid boundaries.
    void setExtrapolator(Extrapolator* epol);

    /// Sets the default extrapolator.
    /// @todo Needed on public API? Just read from grid file.
    void setDefaultExtrapolator();

    /// Check of valid Interpolator.
    bool hasInterpolator() const {
      return _interpolator != 0;
    }

    /// Check of valid Extrapolator.
    bool hasExtrapolator() const {
      return _extrapolator != 0;
    }

    //@}



    /// Loads the given member by path to the member grid file.
    /// @todo Also need loading by filename, and by set name + member ID
    static PDFGrid* load( PDFGrid*, const YAML::Node&, std::ifstream& );



    /// @name Info about the grid, and access to the raw data points
    //@{

    /// Return knot values in x
    const AxisKnots& xKnots() const {
      return _xknots;
    }

    /// Return knot values in Q2
    const AxisKnots& q2Knots() const {
      return _q2knots;
    }

    /// Get the index of the closest x knot row <= x
    size_t xKnotLow(double x) const {
      /// @todo Test for x in grid range
      size_t i = lower_bound(xKnots().begin(), xKnots().end(), x) - xKnots().begin();
      if (i == xKnots().size()-1) --i; // if last row, step back
      return i;
    }

    /// Get the index of the closest Q2 knot column <= q2
    size_t q2KnotLow(double q2) const {
      size_t i = lower_bound(q2Knots().begin(), q2Knots().end(), q2) - q2Knots().begin();
      if (i == q2Knots().size()-1) --i; // if last col, step back
      return i;
    }

    /// Get the raw xf(x,Q2) data points
    const double* rawdata(PID_t id) const {
      if (!hasPID(id)) {
        std::stringstream error;
        error << "Undefined particle ID requested: " << id;
        throw std::runtime_error(error.str());
      }
      return _rawdata.find(id)->second;
    }

    /// @brief Transform a (ix, iQ2) pair into a 1D "raw" index
    size_t rawindex(size_t ix, size_t iq2) const {
      if (ix >= xKnots().size()) throw std::runtime_error("Invalid x index");
      if (iq2 >= q2Knots().size()) throw std::runtime_error("Invalid Q2 index");
      return ix + iq2 * xKnots().size();
    }

    //@}


  private:

    /// Interpolation grid anchor point lists in x and Q2
    AxisKnots _xknots, _q2knots;

    /// Raw data grids, indexed by flavour
    std::map<PID_t, double*> _rawdata;

    /// Associated interpolator
    Interpolator* _interpolator;

    /// Associated extrapolator
    Extrapolator* _extrapolator;

  };


}
