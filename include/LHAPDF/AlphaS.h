#pragma once

#include <LHAPDF/Utils.h>
#include <LHAPDF/Exceptions.h>

namespace LHAPDF {


  /// Contains different schemes for calculating alpha_s
  ///
  /// The design of the AlphaS classes is that they are substitutible
  /// (cf. polymorphism) and are entirely non-dependent on the PDF and Info
  /// objects: hence they can be used by external code that actually doesn't
  /// want to do anything at all with PDFs, but which just wants to do some
  /// alpha_s interpolation.
  class AlphaS {
  public:

    AlphaS() {
      _qmasses[0] = 0.0017;
      _qmasses[1] = 0.0041;
      _qmasses[2] = 0.1;
      _qmasses[3] = 1.29;
      _qmasses[4] = 4.1;
      _qmasses[5] = 172.5;

      _mz = 91.1876;
      _alphas_mz = 0.118;
    }

    double alphaS_Q(double q) { return alphaS_Q2(q*q); }
    virtual double alphaS_Q2(double q2) = 0;

  protected:

    /// Calculate the number of active flavours at energy scale Q2
    int _nf_Q2(double q2);

    /// Calculate the i'th beta function given the number of active flavours
    ///
    /// Currently limited to 0 <= i <= 2.
    double beta(int i, int nf);

    /// Calculate a vector of beta functions given the number of active flavours
    ///
    /// Currently returns a 3-element vector of beta0 -- beta2.
    std::vector<double> betas(int nf);

    /// Masses of quarks in GeV.  Used to calculate the number of quarks that are active at a given energy range Q2
    double _qmasses[6];

    /// Mass of the Z-boson and value of alpha_s(MZ)
    double _mz,_alphas_mz;

    /// LambdaQCD values for 4 and 5 flavour regimes
    double _lambda4, _lambda5;

    /// Order of QCD (expressed as number of loops)
    int _order;

  };



  class AlphaS_Analytic : public AlphaS {
  public:

    double alphaS_Q2(double q2);

  };



  class AlphaS_ODE : public AlphaS {
  public:

    double alphaS_Q2(double q2);

  private:

    /// Calculate the first order derivative, dAlphaS/dQ2, as it appears in differential equation
    double _derivative(double t, double y, const std::vector<double>& beta);

  };



  /// @todo Complete and re-enable
  // /// Interpolate alpha_s from points in Q2
  // class AlphaS_Ipol : public AlphaS {
  // public:

  //   double alphaS_Q2(double q2);

  // };


}
