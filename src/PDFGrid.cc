#include "LHAPDF/PDFGrid.h"
#include "LHAPDF/Interpolator.h"
#include "LHAPDF/factory.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstring>

namespace LHAPDF {


  void PDFGrid::setInterpolator(Interpolator* i) {
    if (_interpolator != 0) {
      if (_allocatedInterpolator) {
        // Delete the allocated interpolator
        /// @todo Why the distinction?
        delete _interpolator;
      } else {
        // Release association with PDFGrid
        _interpolator->bind(*this);
      }

      // Set interpolator
      _interpolator = i;
      _allocatedInterpolator = false;

      // Allow preprocessing of PDFGrid
      _interpolator->bind( *this );
    }


    void PDFGrid::setDefaultInterpolator() {
      /// @todo Make use of setInterpolator to reduce code

      if (_interpolator != 0) {
        if (_allocatedInterpolator) {
          // Delete the allocated interpolator
          delete _interpolator;
        } else {
          // Release association with PDFGrid
          _interpolator->unbind(*this);
        }
      }
      /// @todo Remove noise
      std::cout << "Ejected old interpolator" << std::endl;
      // Use interpolator factory
      /// @todo Need better key name
      const std::string ipolname = set()->getMetadata("SetInterpolator")
        interpolator = createInterpolator(ipolname);
      allocatedInterpolator = true;
      /// @todo Remove noise
      std::cout << "Loaded new interpolator" << std::endl;
    }

    void PDFGrid::setExtrapolator( Extrapolator* e ) {
      if (_extrapolator != 0) {
        if (_allocatedExtrapolator) {
          delete _extrapolator;
        } else {
          _extrapolator->unbind(*this);
        }
      }
      _allocatedExtrapolator = false;
      _extrapolator = e;
      _extrapolator->bind(*this);
    }

    void PDFGrid::setDefaultExtrapolator() {
      if (_extrapolator != 0) {
        if (_allocatedExtrapolator) {
          delete _extrapolator;
        } else {
          extrapolator->unbind(*this);
        }
      }
      // Use Extrapolator factory
      /// @todo Need better key name
      const std::string xpolname = set()->getMetadata("SetExtrapolator")
        _extrapolator = createExtrapolator(xpolname);
      _allocatedExtrapolator = true;
    }

    double PDFGrid::xfxQ2( PID_t id, const X_t x, const Q2_t q2 ) const {
      ///Error checks
      //Physical Ranges
      if (!inPhysicalRangeQ2( x, q2 )) {
        std::stringstream error;
        error << "Unphysical range given: " << x << " " << q2;

        throw std::runtime_error( error.str() );
      }

      //Undefined PIDs
      if (flavors.find( id ) == flavors.end()) {
        std::stringstream error;
        error << "Undefined flavour requested: " << id;

        throw std::runtime_error( error.str() );
      }

      //OUTSIDE IS EXTRAPOLATION : INSIDE IS INTERPOLATION
      if (inRangeQ2( x, q2 )) {
        //Valid Interpolator
        if( !hasInterpolator()) {
          std::stringstream error;
          error << "Undefined interpolator.";

          throw std::runtime_error( error.str() );
        }

        //Invoke interpolator to find value
        return interpolator->interpolateQ2( *this, id, x, q2 );
      }
      else {
        //Valid Extrapolator
        if (!hasExtrapolator()) {
          std::stringstream error;
          error << "Undefined extrapolator.";

          throw std::runtime_error( error.str() );
        }

        //Invoke extrapolator to find value
        return extrapolator->extrapolateQ2( *this, id, x, q2 );
      }
    }

    Idx_t bisectionSearch( const double x, const AxisKnots& knots ) {
      Idx_t low = 0, high = knots.size()-1;

      while ((high - low) > 1) {
        Idx_t middle = low + ( (high - low) >> 1 );

        if (x > knots[middle])
          low = middle;
        else
          high = middle;
      }

      return low;
    }


    /// @todo Use std library lookups (upper_bound, etc.)
    Idx_t binarySearch( const double x, const AxisKnots& points ) {
      Idx_t low = 0, high = points.size()-1;
      while (true) {
        Idx_t middle = low + ((high - low) >> 2);
        if (middle == low || middle == high)
          break;
        if (x > points[middle])
          low = middle;
        else
          high = middle;
      }
      return low;
    }


    void PDFGrid::lookupClosestLow( const X_t x, const Q2_t q2, Idx_t& xidx, Idx_t& q2idx ) const {
      //DO ERROR CHECKS ON X AND Q?
      //TRUST INTERPOLATORS?

      //Bisection search to find closest low limit
      //POSSIBLE: Implement a hunting algorithm to have a guess
      xidx = bisectionSearch( x, xknots );
      q2idx = bisectionSearch( q2, q2knots );
    }


    void PDFGrid::lookupClosest( const X_t x, const Q2_t q2, Idx_t& xidx, Idx_t& q2idx ) const {
      //Binary Search Axes
      xidx = binarySearch( x, xknots );
      q2idx = binarySearch( q2, q2knots );
    }


    inline void tokenize( const std::string& line, AxisKnots& knots ) {
      size_t start = 0;
      size_t pos = 0;
      do {
        pos = line.find_first_of( " ", start );
        knots.push_back( atof( line.substr( start, (pos - start) ).c_str() ) );
        start = pos+1;
      } while (pos != std::string::npos);
    }



    PDFGrid* PDFGrid::load( PDFGrid* grid, const YAML::Node& head, std::ifstream& file ) {
      // Parse grid knots
      for (YAML::Iterator xsit = head["Xs"].begin(); xsit != head["Xs"].end(); ++xsit) {
        double x;
        (*xsit) >> x;
        grid->xknots.push_back(x);
      }

      for (YAML::Iterator q2sit = head["Q2s"].begin(); q2sit != head["Q2s"].end(); ++q2sit) {
        double q2;
        (*q2sit) >> q2;
        grid->q2knots.push_back(q2);
      }

      // Parse grid data
      // Allocate flavor data
      std::vector<PID_t>::const_iterator piditer = grid->set->getFlavours().begin();
      for( ; piditer != grid->set->getFlavours().end(); ++piditer ) {
        grid->flavors[*piditer] = new double[grid->xknots.size()*grid->q2knots.size()];
      }

      // Parse grid lines
      unsigned int cline = 0;
      for(; cline < grid->xknots.size()*grid->q2knots.size(); ++cline ) {
        // Read a line
        if( !file.good() ) {
          std::stringstream error;
          error << "ifstream ran out of data @ " << cline;

          throw std::runtime_error( error.str() );
        }
        std::string line;
        getline( file, line );

        // Parsing individual grid line
        const char* cstr = line.c_str();
        char str[line.size()+1];  /// @todo Not ISO C++
        memcpy( str, cstr, (line.size()+1)*sizeof(char) );

        char *token, *prog;
        unsigned int flavor=0;
        token = strtok_r( str, " ", &prog );
        while( token != NULL ) {
          //Process token
          grid->flavors[grid->set->getFlavours()[flavor]][cline] = atof( token );

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
