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
      /// @todo And delete the data grids!
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
      for (std::map<PID_t, double*>::const_iterator i = _ptdata.begin(); i != _ptdata.end(); ++i) {
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

    /// @brief Set the interpolator by value
    ///
    /// The passed value must be a concrete instantiation of the Interpolator
    /// interface. It will be copied and heap-assigned for use inside this PDFGrid.
    template <typename INTERPOLATOR>
    void PDFGrid::setInterpolator(INTERPOLATOR ipol) {
      _interpolator = new INTERPOLATOR(ipol);
      _interpolator->bind(this);
    }

    /// @brief Set the interpolator by pointer
    ///
    /// This interpolator argument is exactly the one that will be used by this
    /// PDFGrid: it will not be copied and the PDF takes ownership of the
    /// pointer and will delete the interpolator when the PDF goes out of scope.
    ///
    /// @todo Use smart pointers?
    void setInterpolator(Interpolator* ipol) {
      _interpolator = ipol;
      _interpolator->bind(this);
    }

    /// @brief Set the interpolator by name
    ///
    /// Use the interpolator specified by the given name, as passed to the
    /// createInterpolator factory function.
    void setInterpolator(const std::string& ipolname) {
      setInterpolator(createInterpolator(ipolname));
    }

    /// Get the current interpolator (ownership remains with the PDFGrid).
    const Interpolator* interpolator() const {
      return _interpolator;
    }



    /// @brief Set the extrapolator by value
    ///
    /// The passed value must be a concrete instantiation of the Extrapolator
    /// interface. It will be copied and heap-assigned for use inside this PDFGrid.
    template <typename EXTRAPOLATOR>
    void PDFGrid::setExtrapolator(EXTRAPOLATOR xpol) {
      _extrapolator = new EXTRAPOLATOR(xpol);
      _extrapolator->bind(this);
    }

    /// @brief Set the extrapolator by pointer
    ///
    /// This extrapolator argument is exactly the one that will be used by this
    /// PDFGrid: it will not be copied and the PDF takes ownership of the
    /// pointer and will delete the extrapolator when the PDF goes out of scope.
    ///
    /// @todo Use smart pointers?
    void setExtrapolator(Extrapolator* xpol) {
      _extrapolator = xpol;
      _extrapolator->bind(this);
    }

    /// @brief Set the extrapolator by name
    ///
    /// Use the extrapolator specified by the given name, as passed to the
    /// createExtrapolator factory function.
    void setExtrapolator(const std::string& xpolname) {
      setExtrapolator(createExtrapolator(xpolname));
    }

    /// Get the current extrapolator (ownership remains with the PDFGrid).
    const Extrapolator* extrapolator() const {
      return _extrapolator;
    }

    //@}



    /// Loads the given member by path to the member grid file.
    /// @todo Also need loading by filename, and by set name + member ID
    static PDFGrid* load(PDFGrid*, const YAML::Node&, std::ifstream&);



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
    const double* ptdata(PID_t id) const {
      if (!hasPID(id)) {
        std::stringstream error;
        error << "Undefined particle ID requested: " << id;
        throw std::runtime_error(error.str());
      }
      return _ptdata.find(id)->second;
    }

    /// @brief Transform a (ix, iQ2) pair into a 1D "raw" index
    size_t ptindex(size_t ix, size_t iq2) const {
      if (ix >= xKnots().size()) throw std::runtime_error("Invalid x index");
      if (iq2 >= q2Knots().size()) throw std::runtime_error("Invalid Q2 index");
      return ix + iq2 * xKnots().size();
    }

    //@}


  private:

    /// Interpolation grid anchor point lists in x and Q2
    AxisKnots _xknots, _q2knots;

    /// Raw data grids, indexed by flavour
    std::map<PID_t, double*> _ptdata;

    /// Associated interpolator
    Interpolator* _interpolator;

    /// Associated extrapolator
    Extrapolator* _extrapolator;

  };


}
