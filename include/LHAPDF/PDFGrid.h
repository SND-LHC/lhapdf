#pragma once

#include "LHAPDF/Types.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include <vector>
#include <map>
#include <stdlib.h>
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
    /// @todo Shouldn't it be possible to construct without a set pointer?
    PDFGrid(const PDFSet* setp=0)
      : interpolator(0), allocatedInterpolator(false),
        extrapolator(0), allocatedExtrapolator(false)
        _set(setp)
    {  }


    /// Destructor
    PDFGrid::~PDFGrid() {
      // Check if the class created the interpolator
      /// @todo Use shared ptrs instead of this manual allocation stuff
      if (interpolator != 0 && allocatedInterpolator)
        delete interpolator;
      if (extrapolator != 0 && allocatedExtrapolator)
        delete extrapolator;
    }


    /// @brief PDF xf(x,Q2) value
    double xfxQ2( PID_t, double x, double q2 ) const;


    /// @name Information about the interpolation grid
    //@{

    /// Checks whether (X,Q2) is in PDF data.
    ///
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    bool inRangeQ2(double x, double q2) const {
      // Check x is in range
      if (x < xknots[0]) return false;
      if (xknots[xknots.size()-1] < x) return false;
      // Check q2 is in range
      if (q2 < q2knots[0]) return false;
      if (q2knots[q2knots.size()-1] < q2) return false;
      // Else...
      return true;
    }


    /// Look up the lowest knot of the rectangle in which (X,Q2) is located.
    ///
    /// This function is useful for interpolators to be able to look up which patch the (X,Q2) is
    /// in and then interpolating between the returned index and the index+1.
    ///
    /// The index returned will never be on the high edge of the grid knots.
    ///
    /// \param X Momentum fraction coordinate.
    /// \param Q2 Energy scale coordinate.
    /// \param XIdx Index in X axis.
    /// \param Q2Idx Index in Q2 axis.
    ///
    /// @todo Return pair<size_t, size_t>
    /// @todo Need Q2 in name?
    std::pair<size_t, size_t> lookupClosestLow(double x, double q2) const;

    /// Looks up the closest grid knot.
    ///
    /// \param X Momentum fraction
    /// \param Q2 Energy scale
    /// \param XIdx Index in X axis.
    /// \param Q2Idx Index in Q2 axis.
    ///
    /// @todo Symmetry -- how does this differ from ClosestLow and why no ClosestHigh?
    ///
    std::pair<size_t, size_t> lookupClosest(double x, double q2) const;

    /// Check whether the (XIdx,Q2Idx) is valid in this PDF.
    bool isValidIndex(size_t ix, size_t iq2) const;

    /// @brief Transform a (ix, iQ2) pair into a 1D index
    ///
    /// Used to access flavor images at knots in this PDF.
    size_t index(size_t ix, size_t iq2) const {
      // X is along a scanline, Q2 is a scanline
      return xidx + q2idx * xknots.size();
    }

    /// Return knot values in x
    const AxisKnots& getXKnots() const {
      return _xknots;
    }

    /// Return knot values in Q2
    const AxisKnots& getQ2Knots() const {
      return _q2knots;
    }

    /// Sets the Interpolator to be used for interpolating between grid knots.
    void setInterpolator(Interpolator* );

    /// Sets the default interpolator.
    /// @todo Needed on public API? Just read from grid file.
    void setDefaultInterpolator();

    /// Sets the Extrapolator to be used for extrapolating outside the grid boundaries.
    void setExtrapolator( Extrapolator* );

    /// Sets the default extrapolator.
    /// @todo Needed on public API? Just read from grid file.
    void setDefaultExtrapolator();

    /// Check of valid Interpolator.
    bool hasInterpolator() const {
      return (interpolator != 0);
    }

    /// Check of valid Extrapolator.
    bool hasExtrapolator() const {
      return (extrapolator != 0);
    }

    //@}

    /// Get the list of available flavours by PDG ID code.
  vector<PID_t> flavors() const {
    return _flavors;
  }

    /// @todo What is a PIDdata (and rename it...)
    const PIDdata getPIDData(PID_t id) const {
      if (!hasPID(id)) {
        std::stringstream error;
        error << "Undefined PID_t requested: " << id;
        throw std::runtime_error(error.str());
      }
      return flavors().find(id)->second;
    }


    /// Loads the given member by path
    ///
    /// \param Path the file path to the ".LHm" file
    /// \return
    //static PDFGrid* load( const std::string&, const PDFSet& set );
    static PDFGrid* load( PDFGrid*, const YAML::Node&, std::ifstream& );


  private:

    /// Momentum fraction knots
    AxisKnots _xknots;

    /// Squared energy scale knots
    AxisKnots _q2knots;

    /// Flavor map
    PIDmap _flavors;

    /// Associated Interpolator
    Interpolator* _interpolator;
    /// Flag for self allocated Interpolator
    bool _allocatedInterpolator;

    /// Associated Extrapolator
    Extrapolator* _extrapolator;
    /// Flag for self allocated Extrapolator
    bool _allocatedExtrapolator;

  };


}
