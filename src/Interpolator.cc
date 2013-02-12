#include "LHAPDF/Interpolator.h"
#include "LHAPDF/GridPDF.h"

namespace LHAPDF {


double Interpolator::interpolateXQ2(int id, double x, double q2) const {
      // Subgrid lookup
      /// @todo Do this in two stages to cache the KnotArrayNF?
      /// @todo Add flavour error checking
      const KnotArray1F& subgrid = pdf().subgrid(id, q2);
      // Index look-up
      /// @todo Cache this index lookup for performance
      const size_t ix = subgrid.ixbelow(x);
      const size_t iq2 = subgrid.iq2below(q2);
      // Call the overloaded interpolation routine on this subgrid
      return _interpolateXQ2(subgrid, x, ix, q2, iq2);
    }


}
