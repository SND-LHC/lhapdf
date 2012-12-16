#pragma once

#include "LHAPDF/Types.h"
#include "boost/lexical_cast.hpp"
#include <vector>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <map>

namespace LHAPDF {


  // Forward declaration
  class PDFSet;


  /// @brief PDF is the general interface for access to parton density information.
  ///
  /// The PDF interface declares the general form of all PDF types, such as Grid based or analytic.
  class PDF {
  public:

    /// Virtual destructor, to allow inheritance
    virtual ~PDF() { }

    /// Access the PDFSet
    const PDFSet* pdfSet() const {
      return _set;
    }


    /// @name
    //@{

    /// @brief Get the PDF xf(x) value at (x,q2) for the given PID.
    ///
    /// All grids are defined in q2 rather than q since the natural value
    /// in MC programs is squared, so we typically avoid an expensive sqrt() call.
    /// This is the key function to be overridden in concrete PDF types.
    ///
    /// @param id the parton ID in the PDG scheme
    /// @param x the momentum fraction
    /// @param Q2 the energy scale (squared)
    /// @return the value of xf(x,q2)
    ///
    virtual double xfxQ2(PID_t id, double x, double q2) const = 0;


    /// @brief Get the PDF xf(x) value at (x,q) for the given PID.
    ///
    /// xfxQ will square the given q and return the value from xfxQ2.
    /// All grids are defined in q2 rather than q since the natural value
    /// in MC programs is squared, so we typically avoid an expensive sqrt() call.
    ///
    /// @param id the parton ID in the PDG scheme
    /// @param x the momentum fraction
    /// @param Q the energy scale
    /// @return the value of xf(x,Q)
    ///
    double xfxQ(PID_t id, double x, double q) const {
      return xfxQ2(id, x, q*q);
    }


    /// @brief Get the PDF xf(x) value at (x,q2) for all supported PIDs.
    ///
    /// All grids are defined in q2 rather than q since the natural value
    /// in MC programs is squared, so we typically avoid an expensive sqrt() call.
    /// This is the key function to be overridden in concrete PDF types.
    ///
    /// @param x the momentum fraction
    /// @param Q2 the energy scale (squared)
    /// @return the value of xf(x,q2)
    ///
    std::map<PID_t, double> xfxQ2(double x, double q2) const {
      std::map<PID_t, double> rtn;
      foreach (PID_t id, flavors()) {
        rtn[id] = xfxQ2(id, x, q2);
      }
      return rtn;
    }

    /// @brief Get the PDF xf(x) value at (x,q) for all supported PIDs.
    ///
    /// xfxQ will square the given q and return the value from xfxQ2.
    /// All grids are defined in q2 rather than q since the natural value
    /// in MC programs is squared, so we typically avoid an expensive sqrt() call.
    ///
    /// @param x the momentum fraction
    /// @param Q the energy scale
    /// @return the value of xf(x,Q)
    ///
    std::map<PID_t, double> xfxQ(double x, double q) const {
      return xfxQ2(x, q*q);
    }

    //@}


    /// @name Range checks
    //@{

    /// Check whether the given x is physically valid
    ///
    /// Returns false for x less than 0 or greater than 1, since it
    /// is a momentum fraction and not valid outside those values.
    bool inPhysicalRangeX(double x) const {
      return x >= 0.0 && x < 1.0;
    }

    /// Check whether the given Q2 is physically valid
    ///
    /// Returns false for Q2 less than 0 (Q must be real-valued).
    bool inPhysicalRangeQ2(double q2) const {
      return q2 >= 0.0;
    }

    /// Check whether the given Q is physically valid
    ///
    /// Returns false for Q less than 0 (Q must be positive).
    bool inPhysicalRangeQ(double q) const {
      return inPhysicalRangeQ2(q*q);
    }

    /// Check whether the given (x,Q2) is physically valid
    bool inPhysicalRangeXQ2(double x, double q2) const {
      return inPhysicalRangeX(x) && inPhysicalRangeQ2(q2);
    }

    /// Check whether the given (x,Q) is physically valid
    bool inPhysicalRangeXQ(double x, double q) const {
      return inPhysicalRangeX(x) && inPhysicalRangeQ(q);
    }

    /// inRangeQ will return true when given (X,Q) are in the range of this PDF.
    /// It actually squares given Q and returns value from inRangeQ2.
    ///
    /// \param X the momentum fraction
    /// \param Q the energy scale
    /// \return whether q is in range
    virtual bool inRangeQ(double q) const {
      return inRangeQ2(q*q);
    }

    /// Grid range check for Q2
    ///
    /// inRangeQ2 will return true when given Q2 is in the coverage range of this PDF.
    ///
    /// \param Q2 the squared energy scale
    /// \return
    virtual bool inRangeQ2(double q2) const = 0;

    /// Grid range check for x
    ///
    /// inRangex will return true when given x is in the coverage range of this PDF.
    ///
    /// \param x the momentum fraction
    /// \return
    virtual bool inRangeX(double x) const = 0;

    /// Combined range check for x and Q
    virtual bool inRangeXQ(double x, double q) const {
      return inRangeX(x) && inRangeQ(q);
    }

    /// Combined range check for x and Q2
    bool inRangeXQ2(double x, double q2) const {
      return inRangeX(x) && inRangeQ2(q2);
    }

    //@}


    /// @name Parton content
    //@{

    /// Get the list of particle ID codes for flavors in this PDF.
    virtual std::vector<PID_t> flavors() const = 0;

    /// Checks whether @a id is a valid parton for this PDF.
    bool hasFlavor(PID_t id) const {
      std::vector<PID_t> ids = flavors();
      return std::find(ids.begin(), ids.end(), id) != ids.end();
    }

    //@}


    /// @name Metadata
    //@{

    /// Get all metadata as a map
    const std::map<std::string, std::string>& metadata() const {
      return _metadict;
    }

    /// Get all metadata as a map (non-const)
    std::map<std::string, std::string>& metadata() {
      return _metadict;
    }

    /// Get a specific metadata key as a string
    std::string metadata(const std::string& key) const {
      std::map<std::string, std::string>::const_iterator data = _metadict.find(key);
      if (data != _metadict.end()) {
        throw std::runtime_error("Metadata for key: " + key + " not found.");
      }
      return data->second;
    }

    /// Get a specific metadata key, with automatic type conversion
    template <typename T>
    T metadata(const std::string& key) const {
      std::string s = metadata(key);
      return boost::lexical_cast<T>(s);
    }

    /// Get the name of this PDF set member
    std::string name() const {
      return metadata("Name");
    }

    /// Get the ID code of this PDF set member
    size_t memberID() const {
      return metadata<size_t>("ID");
    }

    /// Get the type of PDF (LO, NLO, etc.)
    int qcdOrder() const {
      return metadata<int>("QCDOrder");
    }

    /// Get the type of PDF error set (Hessian, replicas, etc.)
    std::string errorType() const {
      return metadata("ErrorType");
    }

    //@}


    // /// Load PDF
    // /// @todo What?
    // static PDF* load(const PDFSet*, const std::string& );


  protected:

    /// Pointer back to the containing set
    PDFSet* _set;

    /// Metadata dictionary
    std::map<std::string, std::string> _metadict;

  };


}
