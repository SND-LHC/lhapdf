#include "LHAPDF/Info.h"
#include "yaml-cpp/yaml.h"

namespace LHAPDF {


  void Info::load(const path& filepath) {
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


  /// Load properly cascaded info for a PDF member, including fallback to the set info if it exists
  void Info::loadFull(const path& mempath) { //< @todo Need a better method name!
    const path memberdata = findFile(mempath);
    if (memberdata.empty()) throw ReadError("Could not find PDF data file '" + mempath.native() + "'");
    const string memname = memberdata.filename().native(); //< Can use this to alternatively work out the set name...
    const path setdir = memberdata.parent_path();
    const string setname = setdir.filename().native();
    path setinfo = setdir / (setname + ".info");
    if (exists(setinfo)) load(setinfo.native());
    load(memberdata.native()); //< Override set-level info
  }


}
