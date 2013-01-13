#pragma once

#include "LHAPDF/Info.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"
#include "LHAPDF/Exceptions.h"

namespace LHAPDF {


  /// @brief PDF is the general interface for access to parton density information.
  ///
  /// The PDF interface declares the general form of all PDF types, such as Grid based or analytic.
  class PDF {
  // public:
  protected:

    /// Internal convenience typedef for the AlphaS object handle
    // typedef AlphaS* AlphaSPtr;
    typedef auto_ptr<AlphaS> AlphaSPtr;

    /// @name Creation and deletion
    //@{


    /// Default constructor.
    PDF() { }

    /// Constructor from a file path.
    ///
    /// This constructor reads the member, set, and global metadata and is hence
    /// most useful for being called from the constructors of derived PDF types, e.g.
    /// PDFGrid.
    PDF(const std::string& path) {
      _loadInfo(path);
    }

    /// Constructor from a set name and member ID.
    ///
    /// This constructor reads the member, set, and global metadata and is hence
    /// most useful for being called from the constructors of derived PDF types, e.g.
    /// PDFGrid.
    PDF(const std::string& setname, int member) {
      path searchpath = findFile(pdfmempath(setname, member));
      _loadInfo(searchpath.native());
    }

    /// Constructor from an LHAPDF ID code.
    ///
    /// This constructor reads the member, set, and global metadata and is hence
    /// most useful for being called from the constructors of derived PDF types, e.g.
    /// PDFGrid.
    PDF(int lhaid) {
      const pair<string,int> setname_memid = lookupPDF(lhaid);
      if (setname_memid.second == -1)
        throw IndexError("Can't find a PDF with LHAPDF ID = " + to_str(lhaid));
      path searchpath = pdfmempath(setname_memid.first, setname_memid.second);
      _loadInfo(searchpath.native());
    }

    /// Virtual destructor, to allow unfettered inheritance
    virtual ~PDF() { }

    //@}


  protected:

    void _loadInfo(const path& mempath) {
      /// @todo Set the set name member data somewhere near here
      /// @todo Set the set member number, too
      _mempath = findFile(mempath);
      _info.loadFull(mempath);
    }


  public:

    /// @name PDF values
    //@{

    /// @brief Get the PDF xf(x) value at (x,q2) for the given PID.
    ///
    /// All grids are defined in q2 rather than q since the natural value
    /// in MC programs is squared, so we typically avoid an expensive sqrt() call.
    ///
    /// @param id the parton ID in the PDG scheme
    /// @param x the momentum fraction
    /// @param Q2 the energy scale (squared)
    /// @return the value of xf(x,q2)
    ///
    double xfxQ2(int id, double x, double q2) const {
      // Physical x range check
      if (!inPhysicalRangeX(x)) {
        throw RangeError("Unphysical x given: " + to_str(x));
      }
      // Physical Q2 range check
      if (!inPhysicalRangeQ2(q2)) {
        throw RangeError("Unphysical Q2 given: " + to_str(q2));
      }
      // Undefined PIDs
      if (!hasFlavor(id)) {
        throw FlavorError("Undefined flavour requested: " + to_str(id));
      }
      // Call the delegated method in the concrete PDF object to calculate the in-range value
      return _xfxQ2(id, x, q2);
    }


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
    double xfxQ(int id, double x, double q) const {
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
    std::map<int, double> xfxQ2(double x, double q2) const {
      std::map<int, double> rtn;
      foreach (int id, flavors()) {
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
    std::map<int, double> xfxQ(double x, double q) const {
      return xfxQ2(x, q*q);
    }


  protected:

    /// @brief Calculate the PDF xf(x) value at (x,q2) for the given PID.
    ///
    /// This is the key function to be overridden in concrete PDF types, since
    /// it actually does the calculation of xf(x,Q2) by analytic, interpolation,
    /// or other means. The user-called xfxQ2 method exists so that the physical
    /// range and PID checks need only be done in one place, rather than need to
    /// be re-implemented in each concrete implementation.
    ///
    /// @param id the parton ID in the PDG scheme
    /// @param x the momentum fraction
    /// @param Q2 the energy scale (squared)
    /// @return the value of xf(x,q2)
    ///
    virtual double _xfxQ2(int id, double x, double q2) const = 0;

    //@}


  public:

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


    /// @name Generic metadata
    //@{

    /// Get the info class that actually stores and handles the metadata
    Info& info() { return _info; }

    /// Get the info class that actually stores and handles the metadata (const version)
    const Info& info() const { return _info; }

    /// @todo Set name should be set by the constructor, not the metadata system
    // /// PDF set name
    // std::string setname() const {
    //   return _setname;
    // }

    /// @todo Member ID should be set by the constructor, not the metadata system
    // /// PDF member ID number
    // int memberid() const {
    //   return _memberid;
    // }

    /// Description of the set
    std::string description() const {
      return info().metadata("SetDesc");
    }

    /// Get the type of PDF error set (Hessian, replicas, etc.)
    std::string errorType() const {
      return to_lower_copy(info().metadata("ErrorType"));
    }

    //@}


    /// @name Parton content and QCD parameters
    //@{

    /// List of flavours defined by this PDF set.
    ///
    /// This list is stored locally after its initial read from the Info object
    /// to avoid unnecessary lookups and string decoding, since e.g. it is
    /// looked at by every call to the PDFGrid's Interpolator and Extrapolator
    /// classes.
    const std::vector<int>& flavors() const {
      if (_flavors.empty())
        _flavors = info().metadata< vector<int> >("Flavors");
      return _flavors;
    }

    /// Checks whether @a id is a valid parton for this PDF.
    virtual bool hasFlavor(int id) const {
      const vector<int>& ids = flavors();
      return find(ids.begin(), ids.end(), id) != ids.end();
    }

    /// Order of QCD at which this PDF has been constructed
    ///
    /// "Order" is defined here and throughout LHAPDF as the maximum number of
    /// loops included in the matrix elements, in order to have an integer value
    /// for easy use in comparisons, as opposed to "LO", "NLO", etc. strings.
    int qcdOrder() const {
      return info().metadata<int>("QcdOrder");
    }

    /// @brief Value of alpha_s(Q2) used by this PDF set.
    ///
    /// Calculated numerically, analytically, or interpolated according to
    /// metadata, using the AlphaS classes.
    double alphaS(double q2) const {
      if (_alphas.get() == 0) {
        AlphaS* as = mkAlphaS(info().metadata("AlphaS_Type"));
        /// @todo Configure the QCD params on this AlphaS
        /// @todo Throw an error if the QCD params are changed after a first alpha_s query? How?
        _alphas.reset(as);
      }
      return _alphas->alphaS_Q2(q2);
    }

    //@}


  protected:

    /// Member data file path
    path _mempath;

    /// Metadata container
    Info _info;

    /// Locally cached list of supported PIDs
    mutable vector<int> _flavors;

    /// Optionally loaded AlphaS object
    mutable AlphaSPtr _alphas;

  };


}
