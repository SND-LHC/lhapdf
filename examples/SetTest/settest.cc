#include "PDFSet.h"
#include <iostream>

using namespace LHAPDF;

int main(int argc, const char * argv[]) {
	PDFSet* set = PDFSet::loadByName( "MSTW2008lo90cl_nf3" );
	PDF* mem0 = set->getMember( 0 );
	
	delete set;
}