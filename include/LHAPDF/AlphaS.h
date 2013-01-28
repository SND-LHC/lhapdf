#pragma once

#include "LHAPDF/Utils.h"
#include "LHAPDF/Exceptions.h"

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

    /// Base class constructor for default param setup
    AlphaS();

    /// Calculate alphaS(Q)
    double alphasQ(double q) const { return alphasQ2(q*q); }

    /// Calculate alphaS(Q2)
    /// @todo Throw error in this base method if Q < Lambda?
    virtual double alphasQ2(double q2) const = 0;

    /// Calculate the number of active flavours at energy scale Q
    int nf_Q(double q) const { return nf_Q2(q*q); }

    /// Calculate the number of active flavours at energy scale Q2
    int nf_Q2(double q2) const;

    /// Calculate the i'th beta function given the number of active flavours
    ///
    /// Currently limited to 0 <= i <= 2.
    double beta(int i, int nf) const;

    /// Calculate a vector of beta functions given the number of active flavours
    ///
    /// Currently returns a 3-element vector of beta0 -- beta2.
    std::vector<double> betas(int nf) const;

    /// Get quark masses by PDG code
    double qmass(int id) const { return _qmasses[abs(id)-1]; }
    /// Set quark masses by PDG code
    void setQmass(int id, double value) { _qmasses[abs(id)-1] = value; }

  public:

    /// @name Public properties -- no get/set methods needed
    //@{

    /// Mass of the Z-boson in GeV
    double mz;
    /// Value of alpha_s(MZ)
    double alphas_mz;

    /// LambdaQCD value for 4 flavour regime
    double lambda4;
    /// LambdaQCD value for 5 flavour regime
    double lambda5;

    /// Order of QCD (expressed as number of loops)
    int qcdorder;

    //@}

  protected:

    /// Masses of quarks in GeV.  Used to calculate the number of quarks that are active at a given energy range Q2
    double _qmasses[6];

  };



  /// Calculate alpha_s(Q2) by an analytic approximation
  class AlphaS_Analytic : public AlphaS {
  public:
    double alphasQ2(double q2) const;
  };


  /// Solve the differential equation in alphaS using an implementation of RK4
  class AlphaS_ODE : public AlphaS {
  public:
    double alphasQ2(double q2) const;
  };


  /// Interpolate alpha_s from tabulated points in Q2 via metadata
  class AlphaS_Ipol : public AlphaS {
  public:
    double alphasQ2(double q2) const;
  };


}
