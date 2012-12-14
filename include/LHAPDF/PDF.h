#pragma once

/// @todo Remove/replace?
#include "LHAPDF/types.h"
#include <vector>
#include <sstream>
#include <stdexcept>
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
    const PDFSet* pdfSet() {
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
    std::map<PID_t, double> double xfxQ2(double x, double q2) const {
      std::map<PID_t, double> rtn;
      foreach (PID_t id, pids()) {
        rtn[id] = xfxQ2(id, x, q2);
      }
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

    /// Checks whether (X,Q2) is a valid physical point.
    ///
    /// This will return false in the following cases:
    ///  1) on points below 0.0 in X, since a negative momentum fraction is not valid.
    ///  2) on points above 1.0 in X, since X is a fraction of momentum.
    ///  3) on points below 0.0 in Q2, TODO
    ///
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    ///
    /// @todo Move the physical range check into LHAPDF::PDF
    /// @todo Better name
    bool inPhysicalRangeQ2(double x, double q2) const {
      if (x < 0.0 || x > 1.0) return false;
      return (q2 >= 0.0);
    }

    /// inRangeQ will return true when given (X,Q) are in the range of this PDF.
    /// It actually squares given Q and returns value from inRangeQ2.
    ///
    /// \param X the momentum fraction
    /// \param Q the energy scale
    /// \return whether q is in range
    virtual inline bool inRangeQ(double q) const {
      return inRangeQ2( x, q*q );
    }

    /// Range check for (X,Q2)
    ///
    /// inRangeQ2 will return true when given (X,Q2) are in the range of this PDF.
    ///
    /// \param X the momentum fraction
    /// \param Q2 the squared energy scale
    /// \return
    virtual bool inRangeQ2(double q2) const {
      return true;
    }

    virtual bool inRangeX(double x) const {
      return true;
    }

    virtual bool inRangeXQ(double x, double q) const {
      return inRangeX(x) && inRangeQ(q);
    }

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
      return flavors().find(id) != flavors().end()
    }

    //@}


    /// @name Metadata
    //@{

    /// Get all metadata as a map
    const map<std::string, std::string>& metadata() const {
      return _metadict;
    }

    /// Get all metadata as a map (non-const)
    map<std::string, std::string>& metadata() {
      return _metadict;
    }

    /// Get a specific metadata key as a string
    std::string metadata(const std::string& key) const {
      std::map<std::string, std::string>::const_iterator data = metadict.find(key);
      if (data != metadict.end()) {
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
    std::string PDF::name() const {
      return metadata("Name");
    }

    /// Get the ID code of this PDF set member
    Member_t PDF::memberID() const {
      return metadata<Member_t>("ID");
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


    /// Load PDF
    /// @todo What?
    static PDF* load(const PDFSet*, const std::string& );


  protected:

    /// Pointer back to the containing set
    const PDFSet* _set;

    /// Metadata dictionary
    std::map<std::string, std::string> _metadict;

  };


}
