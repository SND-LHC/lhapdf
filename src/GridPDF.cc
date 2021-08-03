// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/FileIO.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstring>

using namespace std;

namespace LHAPDF {


  void GridPDF::setInterpolator(Interpolator* ipol) {
    _interpolator.reset(ipol);
    _interpolator->bind(this);
  }

  void GridPDF::setInterpolator(const std::string& ipolname) {
    setInterpolator(mkInterpolator(ipolname));
  }

  void GridPDF::_loadInterpolator() {
    const string ipolname = info().get_entry("Interpolator");
    /// @todo What if there is no Interpolator key?
    setInterpolator(ipolname);
  }

  const Interpolator& GridPDF::interpolator() const {
    if (!hasInterpolator()) throw Exception("No Interpolator pointer set");
    return *_interpolator;
  }



  void GridPDF::setExtrapolator(Extrapolator* xpol) {
    _extrapolator.reset(xpol);
    _extrapolator->bind(this);
  }

  void GridPDF::setExtrapolator(const std::string& xpolname) {
    setExtrapolator(mkExtrapolator(xpolname));
  }

  void GridPDF::_loadExtrapolator() {
    const string xpolname = info().get_entry("Extrapolator");
    /// @todo What if there is no Extrapolator key?
    setExtrapolator(xpolname);
  }

  const Extrapolator& GridPDF::extrapolator() const {
    if (!hasExtrapolator()) throw Exception("No Extrapolator pointer set");
    return *_extrapolator;
  }

  /*
  const KnotArrayNF& GridPDF::subgrid(double q2) const {
    assert(q2 >= 0);
    assert(!q2Knots().empty());
    map<double, KnotArrayNF>::const_iterator it = _knotarrays.upper_bound(q2);
    if (it == _knotarrays.begin())
      throw GridError("Requested Q2 " + to_str(q2) + " is lower than any available Q2 subgrid (lowest Q2 = " + to_str(q2Knots().front()) + ")");
    if (it == _knotarrays.end() && q2 > q2Knots().back())
      throw GridError("Requested Q2 " + to_str(q2) + " is higher than any available Q2 subgrid (highest Q2 = " + to_str(q2Knots().back()) + ")");
    --it; // upper_bound (and lower_bound) returns the entry *above* q2: we need to decrement by one element
    // std::cout << "Using subgrid #" << std::distance(_knotarrays.begin(), it) << std::endl;
    return it->second;
  }
  */

  // MK: translate?
  /*
  const vector<double>& GridPDF::q2Knots() const {
    if (_q2knots.empty()) {
      // Get the list of Q2 knots by combining all subgrids
      for (const auto& q2_ka : _knotarrays) { //< auto to make clang happy re. key constness
        const KnotArrayNF& subgrid = q2_ka.second;
        const KnotArray1F& grid1 = subgrid.get_first();
        if (grid1.q2s().empty()) continue; //< @todo This shouldn't be possible, right? Throw instead, or ditch the check?
        for (double q2 : grid1.q2s()) {
          if (_q2knots.empty() || q2 != _q2knots.back()) _q2knots.push_back(q2);
        }
      }
    }
    return _q2knots;
  }
  */


  double GridPDF::_xfxQ2(int id, double x, double q2) const {
    /// Decide whether to use interpolation or extrapolation... the sanity checks
    /// are done in the public PDF::xfxQ2 function.
    //cout << "From GridPDF[0]: x = " << x << ", Q2 = " << q2 << endl;
    double xfx = 0;
    // MK: write propper function
    //int _id = data._pidLookup.find(id)->second;
    int _id;
    if(id != 21){
      _id = data._lookup[id + 6];
    } else{
      _id = data._lookup[0 + 6];
    }
    if(_id == -1) return 0;
    
    if (inRangeXQ2(x, q2)) {
      // cout << "From GridPDF[ipol]: x = " << x << ", Q2 = " << q2 << endl;
      // cout << "Num subgrids = " << _knotarrays.size() << endl;
      // int i = 0;
      // for (std::map<double, KnotArrayNF>::const_iterator it = _knotarrays.begin(); it != _knotarrays.end(); ++it)
      //   cout << "#" << i++ << " from Q = " << sqrt(it->first) << endl;
      xfx = interpolator().interpolateXQ2(_id, x, q2);
    } else {
      // cout << "From GridPDF[xpol]: x = " << x << ", Q2 = " << q2 << endl;
      xfx = extrapolator().extrapolateXQ2(_id, x, q2);
    }
    return xfx;
  }
  
  void GridPDF::_xfxQ2(double x, double q2, std::vector<double>& ret) const {
    //std::cout << "GridPDF::_xfxQ2 : " << x << " " << q2 << std::endl;
    if (inRangeXQ2(x, q2)) {
      interpolator().interpolateXQ2(x, q2, ret);
    } else {
      // Originally, this part was done in PDF::xfxQ2(double, double, vector)
      const int _n = flavors().size();
      ret.clear();
      ret.resize(_n - 1);
      for (int i = 0; i < _n; ++i) {
	extrapolator().extrapolateXQ2(_n, x, q2);
      }
    }
  }


  namespace {

    // A wrapper for std::strtod and std::strtol, for fast tokenizing when all
    // input is guaranteed to be numeric (as in this data block). Based very
    // closely on FastIStringStream by Gavin Salam.
    class NumParser {
    public:
      // Constructor from char*
      NumParser(const char* line=0) { reset(line); }
      // Constructor from std::string
      NumParser(const string& line) { reset(line); }

      // Re-init to new line as char*
      void reset(const char* line=0) {
        _next = const_cast<char*>(line);
        _new_next = _next;
        _error = false;
      }
      // Re-init to new line as std::string
      void reset(const string& line) { reset(line.c_str()); }

      // Tokenizing stream operator (forwards to double and int specialisations)
      template<class T> NumParser& operator>>(T& value) {
        _get(value);
        if (_new_next == _next) _error = true; // handy error condition behaviour!
        _next = _new_next;
        return *this;
      }

      // Allow use of operator>> in a while loop
      operator bool() const { return !_error; }

    private:
      void _get(double& x) { x = std::strtod(_next, &_new_next); }
      void _get(float& x) { x = std::strtof(_next, &_new_next); }
      void _get(int& i) { i = std::strtol(_next, &_new_next, 10); } // force base 10!

      char *_next, *_new_next;
      bool _error;
    };

  }


  void GridPDF::_loadData(const std::string& mempath) {
    string line, prevline;
    int iblock(0), iblockline(0), iline(0), xsize(0);
    vector<double> knots;
    vector<int> pids;
    vector<double> ipid_xfs;

    // MK: do we really need to read the file twice?
    try{
      IFile file(mempath.c_str());
      NumParser nparser; double ftoken; int itoken;
      while (getline(*file, line)) {
        line = trim(line);
	
        // If the line is commented out, increment the line number but not the block line
        iline += 1;
        if (line.find("#") == 0) continue;
        iblockline += 1;

        if (line != "---") { // if we are not on a block separator line...
          // Block 0 is the metadata, which we ignore here
          if (iblock == 0) continue;
	  
          // Parse the data lines
          nparser.reset(line);
          if (iblockline == 1) { // x knots line
	    if(iblock == 1){
	      while (nparser >> ftoken) knots.push_back(ftoken);
	      if (knots.empty())
		throw ReadError("Empty x knot array on line " + to_str(iline));
	      xsize = knots.size();
	    } else { // the x grid should be the same as for the fist i block
	      int tmp = 0;
	      while (nparser >> ftoken) {
		if(ftoken != knots[tmp])
		  throw ReadError("Mismatch in the x-knots");
		++tmp;
	      }
	    }
	    
          } else if (iblockline == 2) { // Q knots line
            while (nparser >> ftoken) knots.push_back(ftoken*ftoken); // note Q -> Q2
            if (knots.size() == xsize)
              throw ReadError("Empty Q knot array on line " + to_str(iline));
          } else if (iblockline == 3) { // internal flavor IDs ordering line
	    if(iblock == 1){
	      while (nparser >> itoken) pids.push_back(itoken);
	    } else {
	      int tmp = 0;
	      while (nparser >> itoken) {
		if (itoken != pids[tmp])
		  throw ReadError("Mismatch in the pids");
		++tmp;
	      }
	    }
            // Check that each line has many tokens as there should be flavours
            if (pids.size() != flavors().size())
              throw ReadError("PDF grid data error on line " + to_str(iline) + ": " + to_str(pids.size()) +
                              " parton flavors declared but " + to_str(flavors().size()) + " expected from Flavors metadata");
            /// @todo Handle sea/valence representations via internal pseudo-PIDs
	    //  MK: What?
          }
	} else{
	  ++iblock;
	  iblockline = 0;
	}
      }
    } catch (Exception& e) {
      throw;
    } catch (std::exception& e) {
      throw ReadError("Read error while parsing " + mempath + " as a GridPDF data file");
    }

    iblock = 0; iblockline = 0; iline = 0;

    // feed data into KnotArray
    // MK: write proper setter functions
    data._knots = knots;
    for(double knot : knots){
      data._log_knots.push_back(log(knot));
    }
    
    data.shape.resize(3);
    data.shape[0] = xsize;
    data.shape[1] = knots.size() - xsize;
    data.shape[2] = pids.size();
    data._pids = pids;

    // create lookuptable to get index id from pid
    data.initPidLookup();
    
    // sets size of data vector
    ipid_xfs.resize(data.shape[0] * data.shape[1] * data.shape[2]);
    
    int qloc(0), qtot(0);
    try {
      int index(0);
      int xindex(0);
      
      IFile file(mempath.c_str());
      NumParser nparser; double ftoken; int itoken;
      while (getline(*file, line)) {

        // Trim the current line to ensure that there is no effect of leading spaces, etc.
        line = trim(line);
        prevline = line; // used to test the last line after the while loop fails

        // If the line is commented out, increment the line number but not the block line
        iline += 1;
        if (line.find("#") == 0) continue;
        iblockline += 1;

        if (line != "---") { // if we are not on a block separator line...
          // Block 0 is the metadata, which we ignore here
          if (iblock == 0) continue;
          nparser.reset(line);
	  if (iblockline == 2) { // Find out how many q values are there
	    qloc = 0;
            while (nparser >> ftoken) ++qloc;
	  } else if (iblockline < 4){
	    continue;
          } else {
            while (nparser >> ftoken) {
	      ipid_xfs[xindex*data.shape[2]*data.shape[1]
		       + qtot*data.shape[2]
		       + index] = ftoken;
	      ++index;
            }	    
	    if( (iblockline != 3) && (iblockline - 3) % qloc == 0){
	      ++xindex;
	      index = 0;
	    }
            // Check that each line has many tokens as there should be flavours
            if (index % pids.size() != 0)
	      // MK: Error message gives wrong output bc. index % pids.size() is not the number of pids
              throw ReadError("PDF grid data error on line " + to_str(iline) + ": " + to_str(index % pids.size()) +
                              " flavor entries seen but " + to_str(pids.size()) + " expected");
          }

        } else { // we *are* on a block separator line	
          // Check that the expected number of data lines were seen in the last block
	  // MK: Does not work anymore, how to translate?
	  /*
          if (iblock > 0 && iblockline - 1 != int(xs.size()*q2s.size()) + 3)
            throw ReadError("PDF grid data error on line " + to_str(iline) + ": " +
                            to_str(iblockline-1) + " data lines were seen in block " + to_str(iblock-1) +
                            " but " + to_str(xs.size()*q2s.size() + 3) + " expected");
	  */
			    
          // Ignore block registration if we've just finished reading the 0th (metadata) block
          if (iblock > 0) {
	    
            // Throw if the last subgrid block was of zero size
            if (ipid_xfs.empty())
              throw ReadError("Empty xf values array in data block " + to_str(iblock) + ", ending on line " + to_str(iline));
          }

          // Increment/reset the block and line counters, etc
          iblock += 1;
          iblockline = 0;
	  index = 0;
	  xindex = 0;
	  qtot += qloc;
        }
	
      }
      // MK: write proper setter methods 
      data._grid = ipid_xfs;
            
      // File reading finished: complain if it was not properly terminated
      if (prevline != "---")
        throw ReadError("Grid file " + mempath + " is not properly terminated: .dat files MUST end with a --- separator line");
      // Error handling
    } catch (Exception& e) {
      throw;
    } catch (std::exception& e) {
      throw ReadError("Read error while parsing " + mempath + " as a GridPDF data file");
    }

    // precompute derivative values
    // MK: use proper setter/getter
    data._dgrid.resize(data.shape[0] * data.shape[1] * data.shape[2]);

    const int nxknots = data.xsize();
    for(int ix(0); ix<nxknots; ++ix){
      for(int iq2(0); iq2<data.q2size(); ++iq2){
	for(int id(0); id<data.size(); ++id){
	  double derivative = 0;
	  // in principle, we could get rid of the if/else's but since this only has to be done
	  //   choose the method that reads better
	  if (ix != 0 && ix != nxknots-1) { //< If central, use the central difference
	    /// @note We evaluate the most likely condition first to help compiler branch prediction
	    const double lddx = (data.xf(ix, iq2, id) - data.xf(ix-1, iq2, id)) / (data.logxs(ix) - data.logxs(ix-1));
	    const double rddx = (data.xf(ix+1, iq2, id) - data.xf(ix, iq2, id)) / (data.logxs(ix+1) - data.logxs(ix));
	    derivative = (lddx + rddx) / 2.0;
	  } else if (ix == 0) { //< If at leftmost edge, use forward difference
	    derivative = (data.xf(ix+1, iq2, id) - data.xf(ix, iq2, id)) / (data.logxs(ix+1) - data.logxs(ix));
	  } else if (ix == nxknots-1) { //< If at rightmost edge, use backward difference
	    derivative = (data.xf(ix, iq2, id) - data.xf(ix-1, iq2, id)) / (data.logxs(ix) - data.logxs(ix-1));
	  } else {
	    throw LogicError("We shouldn't be able to get here!");
	  }
	  data._dgrid[ix*data.shape[1]*data.shape[2] + iq2*data.shape[2] + id] = derivative;
	}
      }
    }
  }


}
