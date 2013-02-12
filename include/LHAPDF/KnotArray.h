#pragma once

#include "LHAPDF/Exceptions.h"
#include "boost/multi_array.hpp"

namespace LHAPDF {


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
      size_t ixbelow(double x) const {
        // Test that x is in the grid range
        if (x < xs().front()) throw GridError("x value " + to_str(x) + " is lower than lowest-x grid point at " + to_str(xs().front()));
        if (x > xs().back()) throw GridError("x value " + to_str(x) + " is higher than highest-x grid point at " + to_str(xs().back()));
        // Find the closest knot below the requested value
        size_t i = upper_bound(xs(), x) - xs().begin();
        return --i; // have to step back to get the knot <= x behaviour
      }


      /// Q2 knot accessor
      const std::vector<double>& q2s() const { return _q2s; }
      /// Q2 knot setter
      /// @todo Reverse the order of lookup here to reverse the order of x and Q2 strides in the data file
      void setq2s(const std::vector<double>& q2s) { _q2s = q2s; _xfs.resize(boost::extents[_xs.size()][_q2s.size()]); }

      /// Get the index of the closest x knot row <= x
      size_t iq2below(double q2) const {
        // Test that x is in the grid range
        if (q2 < q2s().front()) throw GridError("Q2 value " + to_str(q2) + " is lower than lowest-Q2 grid point at " + to_str(q2s().front()));
        if (q2 > q2s().back()) throw GridError("Q2 value " + to_str(q2) + " is higher than highest-Q2 grid point at " + to_str(q2s().back()));
        /// Find the closest knot below the requested value
        size_t i = upper_bound(q2s(), q2) - q2s().begin();
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


}
