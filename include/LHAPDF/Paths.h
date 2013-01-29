#pragma once

#include "LHAPDF/Config.h"
#include "LHAPDF/Utils.h"
#include "boost/filesystem.hpp"

namespace LHAPDF {


  using namespace boost;
  using namespace boost::filesystem;

  /// @name File searching and search path handling functions
  //@{

  /// Get the ordered list of search paths, from $LHAPDF_DATA_PATH and the install location
  std::vector<path> paths();


  /// Set the search paths list as a colon-separated string
  void setPaths(const std::string& pathstr);
  /// Set the search paths list
  inline void setPaths(std::vector<string> paths) {
    setPaths(join(paths, ":"));
  }
  /// Set the search paths list
  inline void setPaths(std::vector<path> paths) {
    vector<string> ss; ss.reserve(paths.size());
    foreach (const path& p, paths) ss.push_back(p.native());
    setPaths(ss);
  }


  /// Prepend to the search paths list
  inline void pathsPrepend(const path& p) {
    vector<path> ps = paths();
    ps.insert(ps.begin(), p);
    ps.pop_back(); //< Discard the auto-added fallback path to the installed data prefix
    setPaths(ps);
  }
  /// Prepend to the search paths list
  inline void pathsPrepend(const std::string& p) {
    pathsPrepend(path(p));
  }


  /// Append to the search paths list
  inline void pathsAppend(const path& p) {
    vector<path> ps = paths();
    ps.pop_back(); //< Discard the auto-added fallback path to the installed data prefix
    ps.push_back(p);
    setPaths(ps);
  }
  /// Append to the search paths list
  inline void pathsAppend(const std::string& p) {
    pathsAppend(path(p));
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

  //@}



  /// @name Functions for handling standard LHAPDF filename structures
  //@{

  inline path pdfmempath(const std::string& setname, int member) {
    const string memname = setname + "_" + to_str_zeropad(member) + ".lha";
    path mempath = setname / memname;
    return findFile(mempath);
  }

  inline path findpdfmempath(const std::string& setname, int member) {
    return findFile(pdfmempath(setname, member));
  }

  inline path pdfsetinfopath(const std::string& setname) {
    const string infoname = setname + ".info";
    path setinfo = setname / infoname;
    return setinfo;
  }

  /// @todo Check that set info and mem=0 file are in same dir?
  inline path findpdfsetinfopath(const std::string& setname) {
    return findFile(pdfsetinfopath(setname));
  }

  //@}


}
