#pragma once

#include "LHAPDF/PDFSet.h"
#include "LHAPDF/BilinearInterpolator.h"
#include "LHAPDF/BicubicInterpolator.h"
#include "LHAPDF/NearestPointExtrapolator.h"
#include <string>
#include <sstream>
#include <stdexcept>

namespace LHAPDF {


  /// Interpolator factory
  static Interpolator* createInterpolator(const std::string& name) {
    /// @todo Convert name to lower case for comparisons
    if (name == "linear")
      return new BilinearInterpolator();
    else if (name == "cubic")
      return new BicubicInterpolator();
    else
      throw std::runtime_error("Undeclared interpolator requested: " + name);
  }


  /// Extrapolator factory
  static Extrapolator* createExtrapolator(const std::string& name) {
    /// @todo Convert name to lower case for comparisons
    if (name == "nearest")
      return new NearestPointExtrapolator();
    else
      throw std::runtime_error("Undeclared extrapolator requested: " + name);
  }


}
