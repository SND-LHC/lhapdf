// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_KnotArray_H
#define LHAPDF_KnotArray_H

#include "LHAPDF/Exceptions.h"
#include "LHAPDF/Utils.h"

namespace LHAPDF {


  /// @brief Internal storage class for PDF data point grids
  ///
  /// We use "array" to refer to the "raw" knot grid, while "grid" means a grid-based PDF.
  /// The "1F" means that this is a single-flavour array
  class KnotArray1F {
  public:

    /// Default constructor just for std::map insertability
    KnotArray1F() {}

    /// Constructor from x and Q2 knot values, and an xf value grid as strided list
    KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots, const std::vector<double>& xfs)
      : _xs(xknots), _q2s(q2knots), _xfs(xfs)
    {
      assert(_xfs.size() == size());
      _syncx();
      _syncq2();
    }

    /// Constructor of a zero-valued array from x and Q2 knot values
    KnotArray1F(const std::vector<double>& xknots, const std::vector<double>& q2knots)
      : _xs(xknots), _q2s(q2knots),
        _xfs(size(), 0.0)
    {
      assert(_xfs.size() == size());
      _syncx();
      _syncq2();
    }


    /// @name x stuff
    ///@{

    /// x knot setter
    ///
    /// @note Also zeros the xfs array, which is invalidated by resetting the x knots
    void setxs(const std::vector<double>& xs) {
      _xs = xs;
      _syncx();
      _xfs = std::vector<double>(size(), 0.0);
    }

    /// Number of x knots
    size_t xsize() const { return _xs.size(); }

    /// x knot accessor
    const std::vector<double>& xs() const { return _xs; }

    /// log(x) knot accessor
    const std::vector<double>& logxs() const { return _logxs; }

    /// Hash comparator
    bool samexs(const KnotArray1F& other) const { return _xgridhash == other._xgridhash; }
    size_t xhash() const { return _xgridhash; }

    /// @brief Get the index of the closest x knot row <= x
    ///
    /// If the value is >= x_max, return i_max-1 (for polynomial spine construction)
    size_t ixbelow(double x) const {
      // Test that x is in the grid range
      if (x < xs().front()) throw GridError("x value " + to_str(x) + " is lower than lowest-x grid point at " + to_str(xs().front()));
      if (x > xs().back()) throw GridError("x value " + to_str(x) + " is higher than highest-x grid point at " + to_str(xs().back()));
      // Find the closest knot below the requested value
      size_t i = upper_bound(xs().begin(), xs().end(), x) - xs().begin();
      if (i == xs().size()) i -= 1; // can't return the last knot index
      i -= 1; // have to step back to get the knot <= x behaviour
      return i;
    }

    ///@}


    /// @name Q2 stuff
    ///@{

    /// Q2 knot setter
    ///
    /// @note Also zeros the xfs array, which is invalidated by resetting the Q2 knots
    void setq2s(const std::vector<double>& q2s) {
      _q2s = q2s;
      _syncq2();
      _xfs = std::vector<double>(size(), 0.0);
    }

    /// Number of Q2 knots
    size_t q2size() const { return _q2s.size(); }

    /// Q2 knot accessor
    const std::vector<double>& q2s() const { return _q2s; }

    /// log(Q2) knot accessor
    const std::vector<double>& logq2s() const { return _logq2s; }

    /// Hash comparator for Q2 knots
    bool sameq2s(const KnotArray1F& other) const { return _q2gridhash == other._q2gridhash; }
    size_t q2hash() const { return _q2gridhash; }

    /// Get the index of the closest Q2 knot row <= q2
    ///
    /// If the value is >= q2_max, return i_max-1 (for polynomial spine construction)
    size_t iq2below(double q2) const {
      // Test that Q2 is in the grid range
      if (q2 < q2s().front()) throw GridError("Q2 value " + to_str(q2) + " is lower than lowest-Q2 grid point at " + to_str(q2s().front()));
      if (q2 > q2s().back()) throw GridError("Q2 value " + to_str(q2) + " is higher than highest-Q2 grid point at " + to_str(q2s().back()));
      /// Find the closest knot below the requested value
      size_t i = upper_bound(q2s().begin(), q2s().end(), q2) - q2s().begin();
      if (i == q2s().size()) i -= 1; // can't return the last knot index
      i -= 1; // have to step back to get the knot <= q2 behaviour
      return i;
    }

    ///@}


    /// @name PDF values at (x, Q2) points
    ///@{

    /// Number of x knots
    size_t size() const { return xsize()*q2size(); }

    /// xf value accessor (const)
    const std::vector<double>& xfs() const { return _xfs; }
    /// xf value accessor (non-const)
    std::vector<double>& xfs() { return _xfs; }
    /// xf value setter
    void setxfs(const std::vector<double>& xfs) { _xfs = xfs; }

    /// Get the xf value at a particular indexed x,Q2 knot
    const double& xf(size_t ix, size_t iq2) const { return _xfs[ix*q2size() + iq2]; }

    ///@}


  private:

    /// Synchronise log(x) array and hash from the x array
    void _syncx() {
      _logxs.resize(_xs.size());
      for (size_t i = 0; i < _xs.size(); ++i) _logxs[i] = log(_xs[i]);
      _xgridhash = _mkhash(_xs);
    }


    /// Synchronise log(x) and log(Q2) arrays from the x and Q2 ones
    void _syncq2() {
      _logq2s.resize(_q2s.size());
      for (size_t i = 0; i < _q2s.size(); ++i) _logq2s[i] = log(_q2s[i]);
      _q2gridhash = _mkhash(_q2s);
    }


    /// Utility function for making a hash code from a vector<double>
    size_t _mkhash(const std::vector<double>& xx) const;

    /// List of x knots
    std::vector<double> _xs;
    /// List of log(x) knots, precomputed for efficiency
    std::vector<double> _logxs;
    /// Hash for this set of x knots
    size_t _xgridhash = 0;

    /// List of Q2 knots
    std::vector<double> _q2s;
    /// List of log(Q2) knots, precomputed for efficiency
    std::vector<double> _logq2s;
    /// Hash for this set of Q2 knots
    size_t _q2gridhash = 0;

    /// List of xf values across the 2D knot array, stored as a strided [ix][iQ2] 1D array
    std::vector<double> _xfs;


  };


  /// @brief A collection of {KnotArray1F}s accessed by PID code
  ///
  /// The "NF" means "> 1 flavour", cf. the KnotArray1F name for a single flavour data array.
  class KnotArrayNF {
  public:

    /// How many {KnotArray1F}s are stored in this container?
    size_t size() const { return _map.size(); }

    /// Is this container empty?
    bool empty() const { return _map.empty(); }

    /// Does this contain a KnotArray1F for PID code @a id?
    bool has_pid(int id) const {
      return _map.find(id) != _map.end();
    }

    /// Get the KnotArray1F for PID code @a id
    const KnotArray1F& get_pid(int id) const {
      if (!has_pid(id)) throw FlavorError("Undefined particle ID requested: " + to_str(id));
      return _map.find(id)->second;
    }

    /// Convenience accessor for any valid subgrid, to get access to the x/Q2/etc. arrays
    const KnotArray1F& get_first() const {
      if (empty()) throw GridError("Tried to access grid indices when no flavour grids were loaded");
      return _map.begin()->second;
    }

    /// Get the KnotArray1F for PID code @a id
    void set_pid(int id, const KnotArray1F& ka) {
      _map[id] = ka;
    }

    /// Indexing operator (non-const)
    KnotArray1F& operator[](int id) { return _map[id]; }

    /// Access the xs array
    const std::vector<double>& xs() const { return get_first().xs(); }
    /// Access the log(x)s array
    const std::vector<double>& logxs() const { return get_first().logxs(); }
    /// Get the index of the closest x knot column <= x (see KnotArray1F)
    size_t ixbelow(double x) const { return get_first().ixbelow(x); }

    /// Access the Q2s array
    const std::vector<double>& q2s() const { return get_first().q2s(); }
    /// Access the log(Q2)s array
    const std::vector<double>& logq2s() const { return get_first().logq2s(); }
    /// Get the index of the closest Q2 knot row <= q2 (see KnotArray1F)
    size_t iq2below(double q2) const { return get_first().iq2below(q2); }

  private:

    /// Storage
    std::map<int, KnotArray1F> _map;

  };
  
  /// @brief Internal storage class for PDF data point grids
  ///
  /// We use "array" to refer to the "raw" knot grid, while "grid" means a grid-based PDF.
  /// The "1F" means that this is a single-flavour array
  class KnotArray{
  public:
    
    /// How many flavours are stored in the grid stored
    size_t size() const { return shape.back(); }

    /// How many x knots are there
    size_t xsize() const { return shape[0]; }

    /// How many q2 knots are there
    size_t q2size() const { return shape[1]; }
    
    /// Is this container empty?
    bool empty() const { return _grid.empty(); }
    
    // General function to find the knot below a given value
    size_t indexbelow(double value, const std::vector<double>& knots) const {
      size_t i = upper_bound(knots.begin(), knots.end(), value) - knots.begin();
      if (i == knots.size()) i -= 1; // can't return the last knot index
      i -= 1;                // step back to get the knot <= x behaviour
      return i;
    }
    
    size_t ixbelow(double x) const {
      return indexbelow(x, _xs);
    }
    
    size_t iq2below(double q2) const {
      return indexbelow(q2, _q2s);
    }
    
    double getxf(int ix, int iq2, int ipid) const {
      return _grid[ix*shape[2]*shape[1] + iq2*shape[2] + ipid];
    }

    double xf(int ix, int iq2, int pid) const {
      return getxf(ix, iq2, pid); 
    }
    
    void setxf(double ix, double iq, int pid, double value){
      _grid[ix*shape[1]*shape[2] + iq*shape[2] + pid] = value;
    }

    double dxf(int ix, int iq2, int ipid) const {
      return _dgrid[ix*shape[2]*shape[1] + iq2*shape[2] + ipid];
    }

    const double& coeff(int ix, int iq2, int pid, int in) const {
      return _coeffs[ix*(shape[1])*shape[2]*4 + iq2*shape[2]*4 + pid*4 + in];
    }
    
    /// Access the Q2s array
    // as they dont exist anymore in the old form, this now returns a vector
    // instead of the reference to one
    const std::vector<double>& xs() const {
      return _xs;
    }

    const std::vector<double>& logxs() const {
      return _logxs;
    }

    
    double xs(int id) const {
      return _xs[id];
    }

    double logxs(int id) const {
      return _logxs[id];
    }
    
    /// Access the Q2s array
    const std::vector<double>& q2s() const {
      return _q2s;
    }
    
    const std::vector<double>& logq2s() const {
      return _logq2s;
    }

    double q2s(int id) const {
      return _q2s[id];
    }

    double logq2s(int id) const {
      return _logq2s[id];
    }
    
    // Access the grid
    const std::vector<double>& grid() const {
      return _grid;
    }

    // Non const access the knots used for the filling
    std::vector<double>& setxknots() {
      return _xs;
    }

    // Non const acess to the q2knots used for the filling
    std::vector<double>& setq2knots() {
      return _q2s;
    }

    std::vector<int>& setPids() {
      return _pids;
    }

    std::vector<size_t>& setShape(){
      return shape;
    }
    
    // Access the polynomial coefficients
    const std::vector<double>& coefficients() const {
      return _coeffs;
    }
    
    bool inRangeX(double x) const {
      if(x < _xs.front()) return false;
      if(x > _xs.back())   return false;
      return true;
    }
    
    bool inRangeQ2(double q2) const {
      if(q2 < _q2s.front()) return false;
      if(q2 > _q2s.back())   return false;
      return true;
    }

    int findPidInPids(int pid) const {
      std::vector<int>::const_iterator it = std::find(_pids.begin(), _pids.end(), pid);
      if(it == _pids.end())
	return -1;
      else
	return std::distance(_pids.begin(), it);
    }

    bool has_pid(int id) const {
      int _id = findPidInPids(id);
      return _id != -1;
    }

    inline int get_pid(int id) const {
      // hardcoded lookup table for particle ids
      // -6,...,-1,21/0,1,...,6,22
      // if id outside of this range, search in list of ids
      if(id < 21) return _lookup[id + 6];
      else if (id == 21) return _lookup[0 + 6];
      else if (id == 22) return _lookup[13];
      else return findPidInPids(id);
    }
    
    void initPidLookup(){
      _lookup.clear();
      if(_pids.size() == 0){
	// MK: is there a better LHAPDF error for that?
	std::cerr << "Internal error when constructing lookup table; need to fill pids before construction"<< std::endl;
	throw;
      }
      for(int i(-6); i<0; i++)
	_lookup.push_back(findPidInPids(i));
      
      _lookup.push_back(findPidInPids(21));
      for(int i(1); i<=6; i++)
	_lookup.push_back(findPidInPids(i));
      _lookup.push_back(findPidInPids(22));
    }

    void fillLogKnots() {
      _logxs.resize(_xs.size());
      for(size_t i(0); i<_xs.size(); ++i)
	_logxs[i] = log(_xs[i]);

      _logq2s.resize(_q2s.size());
      for(size_t i(0); i<_q2s.size(); ++i)
	_logq2s[i] = log(_q2s[i]);
    }
    
    //private:

    // Shape of the interpolation grid
    std::vector<size_t> shape;
    
    // Gridvalues
    std::vector<double> _grid;

    // precompute derivatives
    std::vector<double> _dgrid;

    // Storage for the precomputed polynomial coefficients
    std::vector<double> _coeffs;
    
    // order the pids are filled in
    std::vector<int> _pids;
    std::vector<int> _lookup;

    // knots
    std::vector<double> _xs;
    std::vector<double> _q2s;
    std::vector<double> _logxs;
    std::vector<double> _logq2s;

  };
  

  /// Internal storage class for alpha_s interpolation grids
  class AlphaSArray {
  public:

    /// @name Construction etc.
    ///@{

    /// Default constructor just for std::map insertability
    AlphaSArray() {}

    /// Constructor from Q2 knot values and alpha_s values
    AlphaSArray(const std::vector<double>& q2knots, const std::vector<double>& as)
      : _q2s(q2knots), _as(as)
    {
      _syncq2s();
    }

    ///@}


    /// @name Q2 stuff
    ///@{

    /// Q2 knot vector accessor
    const std::vector<double>& q2s() const { return _q2s; }

    /// log(Q2) knot vector accessor
    const std::vector<double>& logq2s() const { return _logq2s; }

    /// Get the index of the closest Q2 knot row <= q2
    ///
    /// If the value is >= q2_max, return i_max-1 (for polynomial spine construction)
    size_t iq2below(double q2) const {
      // Test that Q2 is in the grid range
      if (q2 < q2s().front()) throw AlphaSError("Q2 value " + to_str(q2) + " is lower than lowest-Q2 grid point at " + to_str(q2s().front()));
      if (q2 > q2s().back()) throw AlphaSError("Q2 value " + to_str(q2) + " is higher than highest-Q2 grid point at " + to_str(q2s().back()));
      /// Find the closest knot below the requested value
      size_t i = upper_bound(q2s().begin(), q2s().end(), q2) - q2s().begin();
      if (i == q2s().size()) i -= 1; // can't return the last knot index
      i -= 1; // have to step back to get the knot <= q2 behaviour
      return i;
    }

    /// Get the index of the closest logQ2 knot row <= logq2
    ///
    /// If the value is >= q2_max, return i_max-1 (for polynomial spine construction)
    size_t ilogq2below(double logq2) const {
      // Test that log(Q2) is in the grid range
      if (logq2 < logq2s().front()) throw GridError("logQ2 value " + to_str(logq2) + " is lower than lowest-logQ2 grid point at " + to_str(logq2s().front()));
      if (logq2 > logq2s().back()) throw GridError("logQ2 value " + to_str(logq2) + " is higher than highest-logQ2 grid point at " + to_str(logq2s().back()));
      /// Find the closest knot below the requested value
      size_t i = upper_bound(logq2s().begin(), logq2s().end(), logq2) - logq2s().begin();
      if (i == logq2s().size()) i -= 1; // can't return the last knot index
      i -= 1; // have to step back to get the knot <= q2 behaviour
      return i;
    }

    ///@}


    /// @name alpha_s values at Q2 points
    ///@{

    /// alpha_s value accessor (const)
    const std::vector<double>& alphas() const { return _as; }
    // /// alpha_s value accessor (non-const)
    // std::vector<double>& alphas() { return _as; }
    // /// alpha_s value setter
    // void setalphas(const valarray& xfs) { _as = as; }

    ///@}


    /// @name alpha_s derivatives vs (log)Q2, useful for interpolation
    ///@{

    /// Forward derivative w.r.t. logQ2
    double ddlogq_forward(size_t i) const {
      return (alphas()[i+1] - alphas()[i]) / (logq2s()[i+1] - logq2s()[i]);
    }

    /// Backward derivative w.r.t. logQ2
    double ddlogq_backward(size_t i) const {
      return (alphas()[i] - alphas()[i-1]) / (logq2s()[i] - logq2s()[i-1]);
    }

    /// Central (avg of forward and backward) derivative w.r.t. logQ2
    double ddlogq_central(size_t i) const {
      return 0.5 * (ddlogq_forward(i) + ddlogq_backward(i));
    }

    ///@}


  private:

    /// Synchronise the log(Q2) array from the Q2 one
    void _syncq2s() {
      _logq2s.resize(_q2s.size());
      for (size_t i = 0; i < _q2s.size(); ++i) _logq2s[i] = log(_q2s[i]);
    }

    /// List of Q2 knots
    std::vector<double> _q2s;
    /// List of log(Q2) knots
    std::vector<double> _logq2s;
    /// List of alpha_s values across the knot array
    std::vector<double> _as;

  };
}
#endif
