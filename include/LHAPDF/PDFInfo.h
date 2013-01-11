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
    void load(const path& filepath);//  {
    //   // Read the YAML file into the metadata map
    //   try {
    //     std::ifstream info(path.c_str());
    //     YAML::Parser parser(info);
    //     YAML::Node doc;
    //     parser.GetNextDocument(doc);
    //     for (YAML::Iterator it = doc.begin(); it != doc.end(); ++it) {
    //       string key, val;
    //       it.first() >> key;
    //       try {
    //         // Assume the value is a scalar type -- it'll throw an exception if not
    //         it.second() >> val;
    //       } catch (const YAML::InvalidScalar& ex) {
    //         // It's a list: process the entries individually into a comma-separated string
    //         string subval;
    //         for (size_t i = 0; i < it.second().size(); ++i) {
    //           it.second()[i] >> subval;
    //           val += subval + ((i < it.second().size()-1) ? "," : "");
    //         }
    //       }
    //       //cout << key << ": " << val << endl;
    //       _metadict[key] = val;
    //     }
    //   } catch (const YAML::ParserException& ex) {
    //     throw ReadError("YAML parse error in " + path + " :" + ex.what());
    //   } catch (const LHAPDF::Exception& ex) {
    //     throw;
    //   } catch (const std::exception& ex) {
    //     throw ReadError("Trouble when reading " + path + " :" + ex.what());
    //   }

    // }


    /// Load properly cascaded info for a PDF member, including fallback to the set info if it exists
    void loadFull(const path& mempath);//  { //< @todo Need a better method name!
    //   const path memberdata = findFile(mempath);
    //   if (memberdata.empty()) throw ReadError("Could not find PDF data file '" + mempath + "'");
    //   const string memname = memberdata.filename().native(); //< Can use this to alternatively work out the set name...
    //   const path setdir = memberdata.parent_path();
    //   const string setname = setdir.filename().native();
    //   path setinfo = setdir / setname + ".info";
    //   if (exists(setinfo)) load(setinfo.native());
    //   load(memberdata.native()); //< Override set-level info
    // }

    //@}


    /// Get the singleton global configuration object
    ///
    /// The global config is populated by reading from lhapdf.conf if it is
    /// found in the search paths.
    static Info& config() {
      static Info _cfg;
      string confpath = findFile("lhapdf.conf").native();
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
    /// Specialisations are defined for unpacking of comma-separated lists of strings, ints, and doubles
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



  // /// Metadata for PDF sets
  /// @todo Re-enable?
  // class SetInfo : public Info {
  // public:
  // };


  // /// Metadata for PDF members
  /// @todo Re-enable?
  // class PDFInfo : public Info {
  // public:
  // };


}


// // Example program
// int main() {
//   using namespace std;
//   const LHAPDF::Info& cfg = LHAPDF::Info::config();
//   cout << cfg.metadata("Foo") << endl;
//   cout << cfg.metadata("Flavors") << endl;
//   vector<int> pids = cfg.metadata< vector<int> >("Flavors");
//   foreach (int f, pids) cout << f << endl;
//   cout << LHAPDF::to_str(pids) << endl;
//   return 0;
// }
