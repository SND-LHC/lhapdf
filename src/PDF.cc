#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"

namespace LHAPDF {


  void PDF::print(std::ostream& os) const {
    os << set().name() << ", member #" << memberID() << ", version " << dataversion() << "\n"
       << "LHAPDF/PDFLIB ID = " << lhapdfID() << "\n"
       << "Flavor content = " << to_str(flavors()) << "\n"
       << description() << std::endl;
  }


}
