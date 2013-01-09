#pragma once

#include "LHAPDF/PDF.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Factories.h"
#include "boost/multi_array.hpp"

namespace LHAPDF {


  // Forward declarations
  class Interpolator;
  class Extrapolator;


  /// @brief A PDF defined via an interpolation grid
  class PDFGrid : public PDF {
  public:

    /// Constructor
    PDFGrid()
      : _interpolator(0), _extrapolator(0)
    {
      /// @todo Parse metadata file, create set if needed, set up alpha_s object, etc.
    }


    /// Destructor
    ~PDFGrid() {
      /// @todo And delete the data grids!
      delete _interpolator;
      delete _extrapolator;
    }


    /// Metadata
    //@{

    /// Get the list of available flavours by PDG ID code.
    /// @todo Or get the flavour list from the set?
    std::vector<int> flavors() const {
      std::vector<int> rtn;
      for (std::map<int, double*>::const_iterator i = _ptdata.begin(); i != _ptdata.end(); ++i) {
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
    void setInterpolator(INTERPOLATOR ipol) {
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
    void setExtrapolator(EXTRAPOLATOR xpol) {
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
    /// @todo Also need loading by set name + member ID
    /// @todo Clarify the ownership of the returned pointer
    static PDFGrid load(const std::string& filepath);
    static PDFGrid load(std::istream& file);


    /// @name Info about the grid, and access to the raw data points
    //@{

    /// Return knot values in x
    const std::vector<double>& xKnots() const {
      return _xknots;
    }

    /// Return knot values in Q2
    const std::vector<double>& q2Knots() const {
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
    const double* ptdata(int id) const {
      if (!hasFlavor(id)) {
        std::stringstream error;
        error << "Undefined particle ID requested: " << id;
        throw FlavorError(error.str());
      }
      return _ptdata.find(id)->second;
    }

    /// @brief Transform a (ix, iQ2) pair into a 1D "raw" index
    size_t ptindex(size_t ix, size_t iq2) const {
      if (ix >= xKnots().size()) throw GridError("Invalid x index");
      if (iq2 >= q2Knots().size()) throw GridError("Invalid Q2 index");
      return ix + iq2 * xKnots().size();
    }

    //@}


  protected:

    /// @brief Get PDF xf(x,Q2) value
    double _xfxQ2(int, double x, double q2) const;


    /// @name Internal storage
    //@{

    /// We use "array" to refer to the "raw" knot grid, while "grid" means a grid-based PDF.
    /// The "1F" means that this is a single-flavour array
    class KnotArray1F {
    public:
      // Use the Boost multi_array for efficiency and ease of indexing
      typedef boost::multi_array<double, 2> valarray;

      KnotArray1F() {} //< for std::map insertability

      KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots, const valarray& xfs)
        : _xs(xknots), _q2s(q2knots), _xfs(xfs)
      { assert(_xfs.shape()[0] == xknots.size() && _xfs.shape()[1] == q2knots.size()); }

      KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots)
        : _xs(xknots), _q2s(q2knots), _xfs(boost::extents[xknots.size()][q2knots.size()])
      { assert(_xfs.shape()[0] == xknots.size() && _xfs.shape()[1] == q2knots.size()); }

      const std::vector<double>& xs() const { return _xs; }
      void setxs(const std::vector<double>& xs) { _xs = xs; _xfs.resize(boost::extents[_xs.size()][_q2s.size()]); }

      const std::vector<double>& q2s() const { return _q2s; }
      void setq2s(const std::vector<double>& q2s) { _q2s = q2s; _xfs.resize(boost::extents[_xs.size()][_q2s.size()]); }

      const valarray& xfs() const { return _xfs; }
      valarray& xfs() { return _xfs; }
      void setxfs(const valarray& xfs) { _xfs = xfs; }

      /// @todo Add index converter and finder methods

    private:
      std::vector<double> _xs, _q2s;
      valarray _xfs;
    };


    /// Typedef for a collection of KnotArray1F accessed by PID code
    /// The "NF" means "> 1 flavour", cf. the KnotArray1F name for a single flavour data array.
    typedef std::map<int, KnotArray1F> KnotArrayNF;

    //@}


  private:

    /// Interpolation grid anchor point lists in x and Q2
    std::vector<double> _xknots, _q2knots;

    /// Raw data grids, indexed by flavour
    /// @todo Need an intermediate type for the subgrids
    std::map<int, double*> _ptdata;

    /// Map of multi-flavour KnotArrays "binned" for lookup by low edge in Q2
    std::map<double, KnotArrayNF> _knotarrays;

    /// Associated interpolator
    Interpolator* _interpolator;

    /// Associated extrapolator
    Extrapolator* _extrapolator;

  };


}
