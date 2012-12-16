#pragma once

#include "LHAPDF/PDF.h"
#include <string>

namespace LHAPDF {


  /// Container of related PDF 'members'
  class PDFSet {
  public:

    /// @todo Need a constructor from a filename

    /// Destructor
    ~PDFSet();

    /// Retrieve general metadata
    std::string metadata( const std::string& key) const {
      std::map<std::string, std::string>::const_iterator data = _metadata.find(key);
      if (data == _metadata.end())
        throw std::runtime_error("Metadata for key: " + key + " not found.");
      return data->second;
    }

    /// Retrieve general metadata, with inline cast
    template <typename T>
    T metadata( const std::string& key) const {
      std::string s = metadata(key);
      return boost::lexical_cast<T>(s);
    }

    /// Set name
    std::string name() const;

    /// Description of the set
    std::string description() const;

    /// Order of QCD at which this PDF has been constructed
    int qcdOrder() const;

    /// @brief Value of alpha_s(Q2) used by this PDF set.
    ///
    /// Calculated numerically, analytically, or interpolated according to metadata.
    double alphaS(double q2) const;

    /// List of flavours defined by this PDF set.
    const std::vector<PID_t>& flavors() const {
      return _flavors;
    }

    /// Get selected member.
    ///
    /// This will load the member if it is not loaded yet or simply return a
    /// reference to the loaded PDF.
    PDF& member(size_t member);

    /// Load the set by path
    static PDFSet* load(const std::string& path);

    /// Load the set by name
    ///
    /// Load a named set in the default search path
    static PDFSet* loadByName( const std::string& );


  private:

    // Declaring unaccessible functions
    /// @todo Why hide these?
    PDFSet(const std::string& path) : _path(path) {}
    PDFSet(const PDFSet& set);
    void operator=(const PDFSet& set);

    /// Holds path for member retrieving
    std::string _path;

    /// Holds all loaded members
    std::map<size_t, PDF*> _members;

    /// Holds all set-level metadata
    std::map<std::string, std::string> _metadata;

    /// Holds all flavours defined by this set
    std::vector<PID_t> _flavors;

  };


}
