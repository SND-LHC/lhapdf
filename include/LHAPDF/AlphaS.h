#pragma once
#include <vector>

namespace LHAPDF {


  /// Contains different schemes for calculating alpha_s
  ///
  /// @todo Rewrite as inheritance from AlphaS base class
  class AlphaS {
  public:
    /// Takes alphaS from metadata
    static double meta(double Q2);

    /// Calculates alphaS using approximate analytical solution to differential equation
    static double analytic(double Q2);

    /// Calculates alphaS using an implementation of RK4 in order to solve differential equation
    /// numerically
    static double numerical(double Q2);


  private:

    /// Calculates the number of active flavours at energy scale Q2
    static int getNfAtQ2(double Q2);

    /// Calculates the first order derivative, dAlphaS/dQ2, as it appears in differential equation
    static double derivative(double t, double y, std::vector<double> beta);

    /// Calculates beta functions given the number of active flavours
    static std::vector<double> getBetas(int nf);

  };


}
