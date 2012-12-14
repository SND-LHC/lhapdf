#pragma once
#include <cstdint>

namespace LHAPDF {


  /// @todo Replace with a string driven by autoconf and move to LHAPDF.h

  static uint8_t version_major = 6;
  static uint8_t version_minor = 0;
  static uint8_t version_rev = 1;

  static uint32_t getVersion( uint32_t size, char* version ) {
    //Check for long enough string
    if ( size < 5)
      return 1;

    //Major
    version[0] = version_major;
    version[1] = '.';

    //Minor
    version[2] = version_minor;
    version[3] = '.';

    //Revision
    version[4] = version_rev;

    return 0;
  }


}
