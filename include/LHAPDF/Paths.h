// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Utils.h"
#include "boost/filesystem.hpp"

namespace LHAPDF {


  using namespace boost;
  using namespace boost::filesystem;

  /// @name File searching and search path handling functions
  //@{

  /// Get the ordered list of search paths, from $LHAPDF_DATA_PATH and the install location
  std::vector<path> paths();
  /// A vector<string> version, for easy Cython compatibility
  inline std::vector<std::string> _paths() {
    std::vector<std::string> rtn;
    foreach (const path& p, paths()) rtn.push_back(p.string());
    return rtn;
  }

  /// Set the search paths list as a colon-separated string
  void setPaths(const std::string& pathstr);
  /// Set the search paths list
  inline void setPaths(std::vector<string> paths) {
    setPaths(join(paths, ":"));
  }
  /// Set the search paths list
  inline void setPaths(std::vector<path> paths) {
    vector<string> ss; ss.reserve(paths.size());
    foreach (const path& p, paths) ss.push_back(p.string());
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
    if (target.empty()) return path();
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
    const string memname = setname + "_" + to_str_zeropad(member) + ".dat";
    path mempath = setname / memname;
    return mempath;
  }

  inline path findpdfmempath(const std::string& setname, int member) {
    return findFile(pdfmempath(setname, member));
  }

  inline path pdfsetinfopath(const std::string& setname) {
    const string infoname = setname + ".info";
    path setinfo = setname / infoname;
    return setinfo;
  }

  inline path findpdfsetinfopath(const std::string& setname) {
    /// @todo Check that set info and mem=0 file are in same dir?
    return findFile(pdfsetinfopath(setname));
  }

  //@}


  /// @brief Get the names of all available PDF sets in the search path
  ///
  /// @note Taken from scanning the directories in the search path
  /// (i.e. LHAPDF_DATA_PATH) for viable PDF sets.
  ///
  /// @note The result is cached when first called, to avoid repeated filesystem
  /// walking. It's assumed that new PDFs will not appear on the filesystem
  /// during a run: please let the authors know if that's not a good assumption!
  inline const std::vector<std::string>& availablePDFSets() {
    static vector<string> rtn;
    if (rtn.empty()) {
      foreach (const path& p, paths()) {
        if (exists(p) && is_directory(p)) {
          directory_iterator idir(p), iend;
          while (idir != iend) {
            const path infopath = idir->path() / (idir->path().filename().string() + ".info");
            if (is_directory(*idir) && exists(infopath)) {
              rtn.push_back(idir->path().filename().string());
            }
            idir++;
          }
        }
      }
      sort(rtn.begin(), rtn.end());
    }
    return rtn;
  }


}
