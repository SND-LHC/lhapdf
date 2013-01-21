#pragma once

#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"
#include "LHAPDF/Exceptions.h"
#include <fstream>

namespace LHAPDF {


  /// Metadata base class for PDFs, PDF sets, or global configuration
  class Info {
  public:

    /// @name Creation and deletion
    //@{

    /// Default constructor
    Info() { }

    /// Constructor
    Info(const std::string& mempath) {
      loadFull(mempath);
    }

    /// Virtual destructor to allow inheritance
    virtual ~Info() { }

    //@}


    /// @name Loading info from YAML files
    //@{

    /// Populate this info object from the specified YAML file path.
    ///
    /// This function may be called several times to read metadata from several
    /// YAML source files. Values for existing keys will be overwritten.
    void load(const path& filepath);

    /// Load properly cascaded info for a PDF member, including fallback to the set info if it exists
    void loadFull(const path& mempath);

    //@}


    /// Get the singleton global configuration object
    ///
    /// @todo Move this to somewhere else -- out of Info?
    ///
    /// The global config is populated by reading from lhapdf.conf if it is
    /// found in the search paths.
    static Info& config() {
      static Info _cfg;
      string confpath = findFile("lhapdf.conf").native();
      // cout << "CONFPATH = " << confpath.empty() << endl;
      if (!confpath.empty()) _cfg.load(confpath);
      return _cfg;
    }


    /// @name General metadata accessors
    //@{

    /// Get all metadata as a map
    const std::map<std::string, std::string>& metadata() const {
      return _metadict;
    }

    /// Get all metadata as a map (non-const)
    std::map<std::string, std::string>& metadata() {
      return _metadict;
    }

    /// Can this Info object return a value for the given key? (it may be defined non-locally)
    bool has_key(const std::string& key) const {
      return has_key_local(key) || Info::config().has_key_local(key);
    }
    /// Is a value defined for the given key on this specific object?
    bool has_key_local(const std::string& key) const {
      return _metadict.find(key) != _metadict.end();
    }

    /// Retrieve a metadata string by key name
    const std::string& metadata(const std::string& key) const {
      if (has_key_local(key)) return _metadict.find(key)->second; //< value is defined locally -- return that
      if (this != &Info::config()) return Info::config().metadata(key); //< if this isn't the global Config, ask that
      throw MetadataError("Metadata for key: " + key + " not found."); //< this is the global Config and key is still not known
    }

    /// Retrieve a metadata string by key name, with an inline type cast
    ///
    /// Specialisations are defined below for unpacking of comma-separated lists
    /// of strings, ints, and doubles.
    template <typename T>
    T metadata(const std::string& key) const {
      const string& s = metadata(key);
      return lexical_cast<T>(s);
    }

    /// Set a keyed value entry
    template <typename T>
    void setMetadata(const std::string& key, const T& val) {
      _metadict[key] = to_str(val);
    }

    //@}


  private:

    /// The string -> string native metadata storage container
    map<string, string> _metadict;

  };


  /// @name Info::metadata function template specialisations
  //@{

  template <>
  inline bool Info::metadata(const std::string& key) const {
    const string& s = metadata(key);
    try {
      bool rtn = lexical_cast<bool>(s);
      return rtn;
    } catch (...) {
      if (s == "true" || s == "on" || s == "yes") return true;
      if (s == "false" || s == "off" || s == "no") return false;
    }
    throw MetadataError("'" + s + "' is not a valid string for conversion to bool type");
  }
  template <>
  inline std::vector<std::string> Info::metadata(const std::string& key) const {
    const string& s = metadata(key);
    vector<string> rtn;
    split(rtn, s, is_any_of(","), token_compress_on);
    return rtn;
  }
  template <>
  inline std::vector<int> Info::metadata(const std::string& key) const {
    const vector<string> strs = metadata< vector<string> >(key);
    vector<int> rtn;
    rtn.reserve(strs.size());
    foreach (const string& s, strs) rtn.push_back( lexical_cast<int>(s) );
    assert(rtn.size() == strs.size());
    return rtn;
  }
  template <>
  inline std::vector<double> Info::metadata(const std::string& key) const {
    const vector<string> strs = metadata< vector<string> >(key);
    vector<double> rtn;
    rtn.reserve(strs.size());
    foreach (const string& s, strs) rtn.push_back( lexical_cast<double>(s) );
    assert(rtn.size() == strs.size());
    return rtn;
  }

  //@}


}
