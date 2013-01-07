#pragma once

#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"
#include "LHAPDF/Exceptions.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

namespace LHAPDF {


  /// Metadata base class for PDFs, PDF sets, or global configuration
  class Info {
  protected:
    Info() { }

  public:

    virtual ~Info() { }

    /// Populate this info object from the specified YAML file path
    void load(const std::string& path) {
      // Read the YAML file into the metadata map
      std::ifstream info(path.c_str());
      YAML::Parser parser(info);
      YAML::Node doc;
      try {
        parser.GetNextDocument(doc);
        for (YAML::Iterator it = doc.begin(); it != doc.end(); ++it) {
          string key, val;
          it.first() >> key;
          try {
            // Assume the value is a scalar type -- it'll throw an exception if not
            it.second() >> val;
          } catch (const YAML::InvalidScalar& ex) {
            // It's a list: process the entries individually into a comma-separated string
            string subval;
            for (size_t i = 0; i < it.second().size(); ++i) {
              it.second()[i] >> subval;
              val += subval + ((i < it.second().size()-1) ? "," : "");
            }
          }
          //cout << key << ": " << val << endl;
          _metadata[key] = val;
        }
      } catch (const YAML::ParserException& ex) {
        cout << "Parse error when reading info from " << path << " :" << ex.what() << endl;
      }
    }


    /// Retrieve general metadata
    const std::string& metadata(const std::string& key) const {
      map<string, string>::const_iterator data = _metadata.find(key);
      if (data == _metadata.end())
        throw MetadataError("Metadata for key: " + key + " not found.");
      return data->second;
    }

    /// Retrieve general metadata, with inline cast
    template <typename T>
    T metadata(const std::string& key) const {
      const string& s = metadata(key);
      return lexical_cast<T>(s);
    }
    template <>
    std::vector<std::string> metadata(const std::string& key) const {
      const string& s = metadata(key);
      vector<string> rtn;
      split(rtn, s, is_any_of(","), token_compress_on);
      return rtn;
    }
    template <>
    std::vector<int> metadata(const std::string& key) const {
      const vector<string> strs = metadata< vector<string> >(key);
      vector<int> rtn;
      rtn.reserve(strs.size());
      foreach (const string& s, strs) rtn.push_back( lexical_cast<int>(s) );
      assert(rtn.size() == strs.size());
      return rtn;
    }
    template <>
    std::vector<double> metadata(const std::string& key) const {
      const vector<string> strs = metadata< vector<string> >(key);
      vector<double> rtn;
      rtn.reserve(strs.size());
      foreach (const string& s, strs) rtn.push_back( lexical_cast<double>(s) );
      assert(rtn.size() == strs.size());
      return rtn;
    }


  private:

    std::map<std::string, std::string> _metadata;

  };


  /// Global LHAPDF config
  class Config : public Info {
  public:

    static Config& make() {
      static Config _cfg;
      /// Read from $prefix/lhapdf.conf if it exists
      string confpath = findFile("lhapdf.conf").native();
      if (!confpath.empty()) _cfg.load(confpath);
      return _cfg;
    }

  };


  /// Metadata for PDF sets
  class SetInfo : public Info {
  public:
  };


  /// Metadata for PDF members
  class PDFInfo : public Info {
  public:

    /// @todo Need a constructor from a filename (and a default one)

    // /// Set name
    // std::string name() const;

    // /// Description of the set
    // std::string description() const;

    // /// Order of QCD at which this PDF has been constructed
    // int qcdOrder() const;

    // /// @brief Value of alpha_s(Q2) used by this PDF set.
    // ///
    // /// Calculated numerically, analytically, or interpolated according to metadata.
    // ///
    // /// @todo Intead return an alpha_s calculator bound to this PDF?
    // double alphaS(double q2) const;

    // /// List of flavours defined by this PDF set.
    // const std::vector<int>& flavors() const;

  };


}


// Example program
int main() {
  using namespace std;
  const LHAPDF::Config& cfg = LHAPDF::Config::make();
  cout << cfg.metadata("Foo") << endl;
  return 0;
}
