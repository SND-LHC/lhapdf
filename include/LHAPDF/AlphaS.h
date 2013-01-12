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

    double alphaS_Q(double q) { return alphaS_Q2(q*q); }
    virtual double alphaS_Q2(double q2) = 0;

  protected:

    /// Calculate the number of active flavours at energy scale Q2
    int _nf_Q2(double q2);

    /// Calculate beta functions given the number of active flavours
    ///
    /// Avoid copying the vector: return by const ref, or just return individual
    /// doubles with an extra index argument?
    std::vector<double> betas(int nf);

    double _mz,_alphas_mz;
    double _lambda4, _lambda5;
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



  /// Interpolate alpha_s from points in Q2
  class AlphaS_Ipol : public AlphaS {
  public:

    double alphaS_Q2(double q2);

  };


}
