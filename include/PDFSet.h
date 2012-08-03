//Include Guard
#pragma once

//LHAPDF Include
#include "types.h"

#include <string>

//Forward Declarations
class PDF;

//PDFSet Declaration
/// The PDFSet is the most general container.

/// A Set contains PDF 'members' and meta data that is either Set general or default meta data for
/// members.
class PDFSet {
public:
    /// Deletes allocated memory when deconstructing.
    ~PDFSet();
    
    /// Get selected member.
    
    /// This will load the member if it is not loaded yet or simply return a pointer to the loaded
    /// instance.
    PDF* getMember( const Member_t );
    
    ///Loads the given member by path
    static PDFSet* load( const std::string& );
    
    ///Loads the given member by name
    static PDFSet* loadByName( const std::string& );

    ///Calculates the renormalised effective QCD coupling constant at Q2 energy scale
    double alphaS(double Q2);
    
private:
    /// Holds all loaded members
    std::map<Member_t, PDF*> members;
};
