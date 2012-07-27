//Include Guard
#pragma once

//Standard Includes
#include <stdint.h>

//Version Declaration & Definition
///
static uint8_t version_major = 6;

///
static uint8_t version_minor = 0;

///
static uint8_t version_rev = 1;

///
static uint32_t getVersion( uint32_t Size, char* Version ) {
    //Check for long enough string
    if( Size < 5)
        return 1;
    
    //Major
    Version[0] = version_major;
    Version[1] = '.';
    
    //Minor
    Version[2] = version_minor;
    Version[3] = '.';
    
    //Revision
    Version[4] = version_rev;

    return 0;
}