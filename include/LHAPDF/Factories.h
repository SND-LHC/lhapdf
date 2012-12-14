#pragma once

#include "LHAPDF/PDFSet.h"
/// @todo need better names with Interpolator/Extrapolator prefixes
#include "LHAPDF/Bilinear.h"
#include "LHAPDF/Bicubic.h"
#include "LHAPDF/Nearest.h"
#include <string>
#include <sstream>
#include <stdexcept>

namespace LHAPDF {


  /// Interpolator factory
  static Interpolator* createInterpolator(const std::string& name) {
    if (name == "bilinear")
      return new BilinearInterpolator();
    else if (name == "bicubic")
      return new BicubicInterpolator();
    else
      throw std::runtime_error("Undeclared interpolator requested: " + name);
  }


  /// Extrapolator factory
  static Extrapolator* createExtrapolator(const std::string& name) {
    if (name == "nearest")
      return new NearestExtrapolator();
    else
      throw std::runtime_error("Undeclared extrapolator requested: " + name);
  }


}
