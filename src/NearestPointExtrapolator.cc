#include "LHAPDF/NearestPointExtrapolator.h"
#include "LHAPDF/GridPDF.h"

namespace LHAPDF {


  namespace { // Unnamed namespace

    // Return the value in the given list that best matches the target value
    double _findClosestMatch(const vector<double>& cands, double target) {
      vector<double>::const_iterator it = lower_bound(cands, target);
      const double upper = *it;
      const double lower = *(--it);
      if (fabs(target - upper) < fabs(target - lower)) return upper;
      return lower;
    }

  }


  double NearestPointExtrapolator::extrapolateXQ2(int id, double x, double q2) const {
    // Find the closest valid x and Q2 points, either on- or off-grid, and use the current interpolator
    const double closestX = (pdf().inRangeX(x)) ? x : _findClosestMatch(pdf().xKnots(), x);
    const double closestQ2 = (pdf().inRangeQ2(q2)) ? q2 : _findClosestMatch(pdf().q2Knots(), q2);
    // cout << "FOO" << closestX << ", " << closestQ2 << endl;
    return pdf().interpolator().interpolateXQ2(id, closestX, closestQ2);
  }


}
