#pragma once

#include <LHAPDF/Config.h>
#include <LHAPDF/Utils.h>
#include <boost/filesystem.hpp>

namespace LHAPDF {


  using namespace boost;
  using namespace boost::filesystem;


  /// Get the ordered list of search paths, from $LHAPDF_DATA_PATH and the install location
  ///
  /// @todo Cache the vector of paths?
  std::vector<path> paths();


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
