#pragma once

#include <LHAPDF/Utils.h>
#include <boost/filesystem.hpp>

namespace LHAPDF {


  using namespace boost;
  using namespace boost::filesystem;


  /// Get the ordered list of search paths, from $LHAPDF_DATA_PATH and the install location
  ///
  /// @todo Cache the vector of paths?
  inline std::vector<path> paths() {
    std::vector<path> rtn;
    // Use LHAPDF_DATA_PATH for all path storage
    char* pathsvar = getenv("LHAPDF_DATA_PATH");
    if (pathsvar != 0) split(rtn, pathsvar, is_any_of(":"), token_compress_on);
    /// @todo Fall back to install prefix
    //pathsvec.push_front(....);
    return rtn;
  }


  /// Return the first location in which a file is found
  ///
  /// If no matching file is found, return an empty path.
  inline path findFile(const path& target) {
    foreach (const path& base, paths()) {
      path p = (!target.is_absolute()) ? base / target : target;
      if (exists(p)) return p;
    }
    return path();
  }


  /// @todo Add functions for pre/appending to and explicitly setting the search paths


}


// Test program:
// #include <iostream>
// int main() {
//   std::cout << 1 << std::endl;
//   foreach (const boost::filesystem::path& p, LHAPDF::paths()) {
//     std::cout << p << std::endl;
//   }
//   std::cout << 2 << std::endl;
//   std::cout << "@" << LHAPDF::findFile("passwd") << "@" << std::endl;
//   return 0;
// }
