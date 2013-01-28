#pragma once

#include "LHAPDF/PDF.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/Extrapolator.h"
#include "boost/multi_array.hpp"

namespace LHAPDF {


  /// @brief A PDF defined via an interpolation grid
  class GridPDF : public PDF {
  public:

    /// @name Creation and deletion
    //@{

    /// Constructor from a file path.
    GridPDF(const std::string& path)
      : PDF(path)
    {
      _loadData(_mempath);
    }

    /// Constructor from a set name and member ID.
    GridPDF(const std::string& setname, int member)
      : PDF(setname, member)
    {
      _loadData(_mempath);
    }

    /// Constructor from a set name and member ID.
    GridPDF(int lhaid)
      : PDF(lhaid)
    {
      _loadData(_mempath);
    }

    /// Virtual destructor to allow inheritance.
    virtual ~GridPDF() { }

    //@}


  protected:

    /// Load the PDF grid data block (not the metadata) from the given PDF member file
    void _loadData(const path& mempath);


  public:

    /// @name Interpolators and extrapolators
    //@{

    /// @brief Set the interpolator by pointer
    ///
    /// The provided Interpolator must have been new'd, as it will not be copied
    /// and ownership passes to this GridPDF: delete will be called on this ptr
    /// when this GridPDF goes out of scope or another setInterpolator call is made.
    void setInterpolator(Interpolator* ipol) {
      _interpolator.reset(ipol);
      _interpolator->bind(this);
    }

    /// @brief Set the interpolator by value
    ///
    /// The passed value must be a concrete instantiation of the Interpolator
    /// interface. It will be copied and heap-assigned for use inside this GridPDF.
    ///
    /// @todo Use SFINAE magic to restrict INTERPOLATOR to subclasses of Interpolator?
    template <typename INTERPOLATOR>
    void setInterpolator(INTERPOLATOR ipol) {
      setInterpolator(new INTERPOLATOR(ipol));
    }

    /// @brief Set the interpolator by name
    ///
    /// Use the interpolator specified by the given name, as passed to the
    /// createInterpolator factory function.
    void setInterpolator(const std::string& ipolname) {
      setInterpolator(mkInterpolator(ipolname));
    }

    /// Find whether an extrapolator has been set on this PDF
    bool hasInterpolator() const {
      return _interpolator.get() != 0;
    }

    /// Get the current interpolator
    const Interpolator& interpolator() const {
      if (_interpolator.get() == 0) { // Load the default interpolator lazily
        // NB. The following is equiv to set-by-name but is explicitly implemented here for const correctness
        const string ipolname = info().metadata("Interpolator");
        Interpolator* ipol = mkInterpolator(ipolname);
        _interpolator.reset(ipol);
        _interpolator->bind(this);
      }
      return *_interpolator;
    }



    /// @brief Set the extrapolator by pointer
    ///
    /// The provided Extrapolator must have been new'd, as it will not be copied
    /// and ownership passes to this GridPDF: delete will be called on this ptr
    /// when this GridPDF goes out of scope or another setExtrapolator call is made.
    void setExtrapolator(Extrapolator* xpol) {
      _extrapolator.reset(xpol);
      _extrapolator->bind(this);
    }

    /// @brief Set the extrapolator by value
    ///
    /// The passed value must be a concrete instantiation of the Extrapolator
    /// interface. It will be copied and heap-assigned for use inside this GridPDF.
    ///
    /// @todo Use SFINAE magic to restrict EXTRAPOLATOR to subclasses of Extrapolator?
    template <typename EXTRAPOLATOR>
    void setExtrapolator(EXTRAPOLATOR xpol) {
      setExtrapolator(new EXTRAPOLATOR(xpol));
    }

    /// @brief Set the extrapolator by name
    ///
    /// Use the extrapolator specified by the given name, as passed to the
    /// createExtrapolator factory function.
    void setExtrapolator(const std::string& xpolname) {
      setExtrapolator(mkExtrapolator(xpolname));
    }

    /// Find whether an extrapolator has been set on this PDF
    bool hasExtrapolator() const {
      return _extrapolator.get() != 0;
    }

    /// Get the current extrapolator
    const Extrapolator& extrapolator() const {
      if (_extrapolator.get() == 0) { // Load the default extrapolator lazily
        // NB. The following is equiv to set-by-name but is explicitly implemented here for const correctness
        const string xpolname = info().metadata("Extrapolator");
        Extrapolator* xpol = mkExtrapolator(xpolname);
        _extrapolator.reset(xpol);
        _extrapolator->bind(this);
      }
      return *_extrapolator;
    }

    //@}


  protected:

    /// @brief Get PDF xf(x,Q2) value (via grid inter/extrapolators)
    double _xfxQ2(int id, double x, double q2) const;


  public:

    /// @name Internal storage
    //@{

    /// @brief Internal storage class for PDF data point grids
    ///
    /// We use "array" to refer to the "raw" knot grid, while "grid" means a grid-based PDF.
    /// The "1F" means that this is a single-flavour array
    class KnotArray1F {
    public:
      /// Use the Boost multi_array for efficiency and ease of indexing
      typedef boost::multi_array<double, 2> valarray;

      /// Default constructor just for std::map insertability
      KnotArray1F() {}

      /// Constructor from x and Q2 knot values, and an xf value grid
      KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots, const valarray& xfs)
        : _xs(xknots), _q2s(q2knots), _xfs(xfs)
      { assert(_xfs.shape()[0] == xknots.size() && _xfs.shape()[1] == q2knots.size()); }

      /// Constructor from x and Q2 knot values
      /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
      KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots)
        : _xs(xknots), _q2s(q2knots), _xfs(boost::extents[xknots.size()][q2knots.size()])
      { assert(_xfs.shape()[0] == xknots.size() && _xfs.shape()[1] == q2knots.size()); }

      /// An explicit operator= is needed due to the Boost multi_array copy semantics
      KnotArray1F& operator=(const KnotArray1F& other) {
        _xs = other._xs;
        _q2s = other._q2s;
        /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
        _xfs.resize(boost::extents[other._xfs.shape()[0]][other._xfs.shape()[1]]);
        _xfs = other._xfs;
        assert(_xfs.shape()[0] == _xs.size() && _xfs.shape()[1] == _q2s.size());
        return *this;
      }

      /// x knot accessor
      const std::vector<double>& xs() const { return _xs; }
      /// x knot setter
      /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
      void setxs(const std::vector<double>& xs) { _xs = xs; _xfs.resize(boost::extents[_xs.size()][_q2s.size()]); }

      /// Get the index of the closest x knot row <= x
      size_t xlow(double x) const {
        // Test that x is in the grid range
        if (x < xs().front()) throw GridError("x value " + to_str(x) + " is lower than lowest-x grid point at " + to_str(xs().front()));
        if (x > xs().back()) throw GridError("x value " + to_str(x) + " is higher than highest-x grid point at " + to_str(xs().back()));
        // Find the closest knot below the requested value
        size_t i = upper_bound(xs().begin(), xs().end(), x) - xs().begin();
        return --i; // have to step back to get the knot <= x behaviour
      }


      /// Q2 knot accessor
      const std::vector<double>& q2s() const { return _q2s; }
      /// Q2 knot setter
      /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
      void setq2s(const std::vector<double>& q2s) { _q2s = q2s; _xfs.resize(boost::extents[_xs.size()][_q2s.size()]); }

      /// Get the index of the closest x knot row <= x
      size_t q2low(double q2) const {
        // Test that x is in the grid range
        if (q2 < q2s().front()) throw GridError("Q2 value " + to_str(q2) + " is lower than lowest-Q2 grid point at " + to_str(q2s().front()));
        if (q2 > q2s().back()) throw GridError("Q2 value " + to_str(q2) + " is higher than highest-Q2 grid point at " + to_str(q2s().back()));
        /// Find the closest knot below the requested value
        size_t i = upper_bound(q2s().begin(), q2s().end(), q2) - q2s().begin();
        return --i; // have to step back to get the knot <= q2 behaviour
      }


      /// xf value accessor (const)
      const valarray& xfs() const { return _xfs; }
      /// xf value accessor (non-const)
      valarray& xfs() { return _xfs; }
      /// xf value setter
      void setxfs(const valarray& xfs) { _xfs = xfs; }

      /// Get the xf value at a particular indexed x,Q2 knot
      /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
      const double& xf(size_t ix, size_t iq2) const { return _xfs[ix][iq2]; }


      // /// @brief Transform a (ix, iQ2) pair into a 1D "raw" index
      // ///
      // /// @todo Rework to use the multi-grid system
      // size_t ptindex(size_t ix, size_t iq2) const {
      //   if (ix >= xKnots().size()) throw GridError("Invalid x index");
      //   if (iq2 >= q2Knots().size()) throw GridError("Invalid Q2 index");
      //   return ix + iq2 * xKnots().size();
      // }


    private:


      /// List of x knots
      vector<double> _xs;
      /// List of Q2 knots
      vector<double> _q2s;
      /// List of xf values across the knot array
      valarray _xfs;

    };


    /// @brief Typedef for a collection of KnotArray1F accessed by PID code
    ///
    /// The "NF" means "> 1 flavour", cf. the KnotArray1F name for a single flavour data array.
    typedef std::map<int, KnotArray1F> KnotArrayNF;

    //@}


    /// @name Info about the grid, and access to the raw data points
    //@{

    /// Get the N-flavour subgrid containing Q2 = q2
    const KnotArrayNF& subgrid(double q2) const {
      assert(q2 >= 0);
      map<double, KnotArrayNF>::const_iterator it = _knotarrays.upper_bound(q2);
      if (it != _knotarrays.begin()) --it; // upper_bound (and lower_bound) return the entry *above* q2: we need to step back
      return it->second;
    }

    /// Get the 1-flavour subgrid for PID=id containing Q2 = q2
    const KnotArray1F& subgrid(int id, double q2) const {
      if (!hasFlavor(id)) throw FlavorError("Undefined particle ID requested: " + to_str(id));
      return subgrid(q2).find(id)->second;
    }

    // /// Return knot values in x
    // const std::vector<double>& xKnots() const {
    //   /// @todo Use Boost join?
    //   return _xknots;
    // }

    // /// Return knot values in Q2
    // const std::vector<double>& q2Knots() const {
    //   /// @todo Use Boost join?
    //   return _q2knots;
    // }


  public:

    /// Check if x is in the grid range
    bool inRangeX(double x) const {
      /// @todo Reinstate this simpler method
      // if (x < xKnots().front()) return false;
      // if (x > xKnots().back()) return false;
      const vector<double>& xs = _knotarrays.begin()->second.begin()->second.xs();
      if (x < xs.front()) return false;
      if (x > xs.back()) return false;
      return true;
    }

    /// Check if q2 is in the grid range
    bool inRangeQ2(double q2) const {
      /// @todo Reinstate this simpler method
      // if (q2 < q2Knots().front()) return false;
      // if (q2 > q2Knots().back()) return false;
      if (q2 < _knotarrays.begin()->first) return false;
      if (q2 > (--_knotarrays.end())->second.begin()->second.q2s().back()) return false;
      return true;
    }

    //@}


  private:

    /// Map of multi-flavour KnotArrays "binned" for lookup by low edge in Q2
    std::map<double, KnotArrayNF> _knotarrays;

    /// Typedef of smart pointer for ipol memory handling
    typedef auto_ptr<Interpolator> InterpolatorPtr;

    /// Typedef of smart pointer for xpol memory handling
    typedef auto_ptr<Extrapolator> ExtrapolatorPtr;

    /// Associated interpolator (mutable to allow laziness)
    mutable InterpolatorPtr _interpolator;

    /// Associated extrapolator (mutable to allow laziness)
    mutable ExtrapolatorPtr _extrapolator;

  };


}
