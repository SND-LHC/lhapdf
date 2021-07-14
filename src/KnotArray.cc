// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/KnotArray.h"
#include <functional>

namespace LHAPDF {


  size_t KnotArray1F::_mkhash(const std::vector<double>& xx) const {
    std::hash<double> hasher;
    size_t rtn = 0;
    for (double x : xx) rtn = 31*rtn + hasher(x);
    return rtn + 1;
  }

  std::vector<size_t> KnotArray::idbelow(std::vector<double> vals){
    std::cerr << "Not implemented yet" << std::endl;
    throw;
  }
  
  const double KnotArray::xf(std::vector<int> ids){
    // Accessor for any size
    // Untested!
    if (ids.size() != shape.size()){
      std::cerr << "Dimensional mismatch" << std::endl;
      throw;
    } else{
      size_t index = 0;
      for(size_t i(0); i<ids.size()-1; ++i){
	size_t blocksize = 1;
	  for(size_t j(i+1); j<ids.size(); ++j)
	    blocksize *= shape[j];
	  index += ids[i]*blocksize;
      }
      index += ids.back();
      return _grid[index];
    }
  }

  

}
