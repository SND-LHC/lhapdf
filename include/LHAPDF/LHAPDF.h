// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

/// @file Just a convenience header that pulls in the main functionality with a single #include

#include "LHAPDF/Version.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Paths.h"

namespace LHAPDF {


  /// @name Global config metadata specialisations
  //@{

  /// @brief All available PDF set names
  ///
  /// @note Taken from scanning the directories in LHAPDF_DATA_PATH for
  /// viable PDF sets.
  const std::vector<std::string>& availablePDFSets() {
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

  //@}


}
