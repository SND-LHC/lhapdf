#include "LHAPDF/Paths.h"
#include "LHAPDF/Info.h"

namespace LHAPDF {


  std::vector<path> paths() {
    vector<path> rtn;
    // Use LHAPDF_DATA_PATH for all path storage
    char* pathsvar = getenv("LHAPDF_DATA_PATH");
    if (pathsvar != 0) split(rtn, pathsvar, is_any_of(":"), token_compress_on);
    // Fall back to current directory if enabled in the central config
    /// @todo Can;t do this, since an infinite loop appears when loading the config()!
    // if (Info::config().has_key("PwdInSearchPath")) {
    //   if (Info::config().metadata<bool>("PwdInSearchPath")) rtn.push_back(".");
    // }
    /// Fall back to install prefix
    rtn.push_back(LHAPDF_DATA_PREFIX);
   return rtn;
  }


  void setPaths(const std::string& pathstr) {
    setenv("LHAPDF_DATA_PATH", pathstr.c_str(), 0);
  }


}
