#pragma once

//Standard Includes
#include <string>
#include <sstream>
#include <exception>

/// UndeclaredInterpolatorException is thrown when the factory tries to instantiate an Interpolator that has not been written.
class UndeclaredInterpolatorException: public std::exception {
public:
    UndeclaredInterpolatorException( const std::string& Name )
    : name( Name ) {}
    
    virtual ~UndeclaredInterpolatorException() throw() {}
    
    const char* what() const throw() {
        std::stringstream message;
        message << "Undeclared interpolator requested of name: " << name << std::endl << ".";
        return message.str().c_str();
    }
    
private:
    const std::string name;
};

/// UndeclaredExtrapolatorException is thrown when the factory tries to instantiate an Extrapolator that has not been written.
class UndeclaredExtrapolatorException: public std::exception {
public:
    UndeclaredExtrapolatorException( const std::string& Name )
    : name( Name ) {}
    
    virtual ~UndeclaredExtrapolatorException() throw() {}
    
    const char* what() const throw() {
        std::stringstream message;
        message << "Undeclared interpolator requested of name: " << name << std::endl << ".";
        return message.str().c_str();
    }
    
private:
    const std::string name;
};

///Include all Interpolator Types
#include "Bilinear.h"
#include "Bicubic.h"

///Interpolator Factory
static Interpolator* createInterpolator( const std::string& Name ) {
    //Insert interpolators here
    if (Name == "bilinear")
        return new BilinearInterpolator();
    else if (Name == "bicubic")
        return new BicubicInterpolator();
    else
        throw UndeclaredInterpolatorException( Name );
}

///Include all Extrapolator Types
#include "Nearest.h"

///Extrapolator Factory
static Extrapolator* createExtrapolator( const std::string& Name ) {
    //Insert extrapolators here
    if (Name == "nearest")
        return new NearestExtrapolator();
    else
        throw UndeclaredExtrapolatorException( Name );
}
