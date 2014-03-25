// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2014 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Paths.h"
#include "LHAPDF/Info.h"

namespace LHAPDF {


  std::vector<std::string> paths() {
    vector<string> rtn;
    // Use LHAPDF_DATA_PATH for all path storage
    char* pathsvar = getenv("LHAPDF_DATA_PATH");
    // But fall back to looking in LHAPATH if the preferred var is not defined
    if (pathsvar == 0) pathsvar = getenv("LHAPATH");
    // Split the paths variable as usual
    if (pathsvar != 0) split(rtn, pathsvar, is_any_of(":"), token_compress_on);
    /// Always look in the install prefix, after other paths are exhausted
    string datadir = LHAPDF_DATA_PREFIX;
    rtn.push_back(datadir / "LHAPDF");
   return rtn;
  }


  void setPaths(const std::string& pathstr) {
    setenv("LHAPDF_DATA_PATH", pathstr.c_str(), 1);
  }


}
