#pragma once

#include "LHAPDF/PDFSet.h"
#include "LHAPDF/BilinearInterpolator.h"
#include "LHAPDF/BicubicInterpolator.h"
#include "LHAPDF/ErrExtrapolator.h"
#include "LHAPDF/NearestPointExtrapolator.h"
#include "boost/algorithm/string.hpp"
#include <string>
#include <sstream>
#include <stdexcept>

namespace LHAPDF {


  /// Interpolator factory
  inline Interpolator* createInterpolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = boost::to_lower_copy(name);
    if (iname == "linear")
      return new BilinearInterpolator();
    else if (iname == "cubic")
      return new BicubicInterpolator();
    else
      throw FactoryError("Undeclared interpolator requested: " + name);
  }


  /// Extrapolator factory
  inline Extrapolator* createExtrapolator(const std::string& name) {
    // Convert name to lower case for comparisons
    const std::string iname = boost::to_lower_copy(name);
    if (iname == "nearest")
      return new NearestPointExtrapolator();
    else if (iname == "error")
      return new ErrExtrapolator();
    else
      throw FactoryError("Undeclared extrapolator requested: " + name);
  }


}
