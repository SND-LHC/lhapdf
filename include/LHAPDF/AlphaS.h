#pragma once

#include <LHAPDF/Utils.h>
#include <LHAPDF/Exceptions.h>

namespace LHAPDF {


  /// Contains different schemes for calculating alpha_s
  ///
  /// @todo Rewrite as inheritance from AlphaS base class
  ///
  /// The design of the AlphaS classes is that they are substitutible
  /// (cf. polymorphism) and are entirely non-dependent on the PDF and PDFInfo
  /// objects: hence they can be used by external code that actually doesn't
  /// want to do anything at all with PDFs, but which just wants to do some
  /// alpha_s interpolation.
  class AlphaS {
  public:

    /// Takes alphaS from metadata
    double meta(double q2);

    /// Calculate alphaS using approximate analytical solution to differential equation
    double analytic(double q2);

    /// Calculate alphaS using an implementation of RK4 in order to solve differential equation
    /// numerically
    double numerical(double q2);

  protected:

    double _lambda4, _lambda5;
    int _order;

  private:

    /// Calculates the number of active flavours at energy scale Q2
    int getNfAtQ2(double q2);

    /// Calculates the first order derivative, dAlphaS/dQ2, as it appears in differential equation
    double derivative(double t, double y, std::vector<double> beta);

    /// Calculates beta functions given the number of active flavours
    std::vector<double> getBetas(int nf);

  };


}
