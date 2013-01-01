#pragma once

#include "LHAPDF/PDF.h"
#include <string>

namespace LHAPDF {


  /// Metadata base class for PDFs, PDF sets, or global configuration
  class Info {
  public:

    /// @todo Need a constructor from a filename (and a default one)

    /// Retrieve general metadata
    std::string metadata( const std::string& key) const {
      std::map<std::string, std::string>::const_iterator data = _metadata.find(key);
      if (data == _metadata.end())
        throw MetadataError("Metadata for key: " + key + " not found.");
      return data->second;
    }

    /// Retrieve general metadata, with inline cast
    template <typename T>
    T metadata( const std::string& key) const {
      std::string s = metadata(key);
      return boost::lexical_cast<T>(s);
    }

    /// @todo Specialise the template on map and list types

  };


  /// Metadata for PDFs
  class PDFInfo : public Info {
  public:

    /// @todo Need a constructor from a filename (and a default one)

    /// Set name
    std::string name() const;

    /// Description of the set
    std::string description() const;

    /// Order of QCD at which this PDF has been constructed
    int qcdOrder() const;

    /// @brief Value of alpha_s(Q2) used by this PDF set.
    ///
    /// Calculated numerically, analytically, or interpolated according to metadata.
    ///
    /// @todo Intead return an alpha_s calculator bound to this PDF?
    double alphaS(double q2) const;

    /// List of flavours defined by this PDF set.
    const std::vector<PID_t>& flavors() const {
      return _flavors;
    }

  };


}
