#include "LHAPDF/Info.h"
#include "LHAPDF/PDFIndex.h"
#include "yaml-cpp/yaml.h"

namespace LHAPDF {


  /// Constructor from a set name and member ID.
  Info::Info(const std::string& setname, int member) {
    path searchpath = findFile(pdfmempath(setname, member));
    loadFull(searchpath.native());
  }

  /// Constructor from an LHAPDF ID code.
  Info::Info(int lhaid) {
    const pair<string,int> setname_memid = lookupPDF(lhaid);
    if (setname_memid.second == -1)
      throw IndexError("Can't find a PDF with LHAPDF ID = " + to_str(lhaid));
    path searchpath = pdfmempath(setname_memid.first, setname_memid.second);
    loadFull(searchpath.native());
  }


  Info& config() {
    static Info _cfg;
    string confpath = findFile("lhapdf.conf").native();
    // cout << "CONFPATH = " << confpath.empty() << endl;
    if (!confpath.empty()) _cfg.load(confpath);
    return _cfg;
  }

  // bool Info::has_key(const std::string& key) const {
  //   return has_key_local(key) || config().has_key_local(key);
  // }

  // const std::string& Info::metadata(const std::string& key) const {
  //   if (has_key_local(key)) return _metadict.find(key)->second; //< value is defined locally -- return that
  //   if (this != &config()) return config().metadata(key); //< if this isn't the global Config, ask that
  //   throw MetadataError("Metadata for key: " + key + " not found."); //< this is the global Config and key is still not known
  // }


  void Info::load(const path& filepath) {
    if (filepath.empty() || !exists(filepath))
      throw ReadError("PDF data file '" + filepath.string() + "' not found");

    // Read the YAML file into the metadata map
    try {
      std::ifstream info(filepath.c_str());
      YAML::Parser parser(info);
      YAML::Node doc;
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
        _metadict[key] = val;
      }
    } catch (const YAML::ParserException& ex) {
      throw ReadError("YAML parse error in " + filepath.native() + " :" + ex.what());
    } catch (const LHAPDF::Exception& ex) {
      throw;
    } catch (const std::exception& ex) {
      throw ReadError("Trouble when reading " + filepath.native() + " :" + ex.what());
    }

  }


  /// @todo Only support loading via PDF set name and member ID, not explicit paths
  /// @todo Replace the loading of the set metadata into the member info with set-level Info singletons
  void Info::loadFull(const path& mempath) { //< @todo Need a better method name!
    // Extract the set name from the member data file path
    const path memberdata = findFile(mempath);
    if (memberdata.empty() || !exists(memberdata)) throw ReadError("Could not find PDF data file '" + mempath.native() + "'");
    const string memname = memberdata.filename().native(); //< Can use this to alternatively work out the set name...
    const path setdir = memberdata.parent_path();
    const string setname = setdir.filename().native();
    path setinfo = findpdfsetinfopath(setname);
    // Load the set info
    if (exists(setinfo)) load(setinfo.native());
    // Load the member info (possibly overriding the set-level metadata)
    load(memberdata.native());
  }


}
