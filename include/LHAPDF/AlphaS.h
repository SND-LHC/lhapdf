// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Utils.h"
#include "LHAPDF/Exceptions.h"

// Riemann zeta function \zeta(3)
// used for calculating the third beta coefficient
#define ZETA_3 1.2020569031

namespace LHAPDF {


  /// @brief Calculator interface for computing alpha_s(Q2) in various ways
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

    /// Destructor
    virtual ~AlphaS() {};

    /// Calculate alphaS(Q)
    double alphasQ(double q) const { return alphasQ2(q*q); }

    /// Calculate alphaS(Q2)
    /// @todo Throw error in this base method if Q < Lambda?
    virtual double alphasQ2(double q2) const = 0;

    /// Calculate the number of active flavours at energy scale Q
    int nf_Q(double q) const { return nf_Q2(q*q); }

    /// Calculate the number of active flavours at energy scale Q2
    virtual int nf_Q2(double q2) const;

    /// Set quark masses by PDG code
    void setQmass(int id, double value);

    /// Set a quark mass (without PDG code).
    /// This allows for more/less than six quarks,
    /// but can not be used together with the other
    /// way of setting quark masses
    void setQmass(double value);

    /// Get the implementation type of this AlphaS
    virtual std::string type() const = 0;

    /// Set the order of QCD (expressed as number of loops)
    void setQCDorder(int order) { _qcdorder = order - 1; }

    // These are defined to avoid problems with
    // for example mkAlphaS... Is there a nicer way to do this?
    // And if not, should these throw an exception/compile time error
    // if used when they are not needed instead of just doing nothing?
    virtual void setMZ(double mz) {};
    virtual void setAlphaSMZ(double alphas) {};
    virtual void setLambda(unsigned int i, double lambda) {};

    //added for AlphaS_Ipol
    virtual void setQ2Values( std::vector<double> ) {};
    virtual void setAlphaSValues( std::vector<double> ) {};

  protected:

    /// Masses of quarks in GeV.  Used to calculate the number
    /// of quarks that are active at a given energy range Q2
    vector<double> _qmasses;
    /// This makes sure you can't mix the two ways of setting masses
    int _masstype;
    /// Calculate the i'th beta function given the number of active flavours
    /// Currently limited to 0 <= i <= 3
    /// Calculated using the MSbar scheme
    double _beta(int i, int nf) const;
    /// Calculate a vector of beta functions given the number of active flavours
    /// Currently returns a 3-element vector of beta0 -- beta3
    std::vector<double> _betas(int nf) const;
    /// Get quark masses by PDG code
    double _qmass(int id) const { return _qmasses[abs(id)-1]; }
    /// Order of QCD (expressed as number of loops)
    int _qcdorder;
  };



  /// Calculate alpha_s(Q2) by an analytic approximation
  class AlphaS_Analytic : public AlphaS {
  public:
    std::string type() const { return "analytic"; }
    double alphasQ2(double q2) const;
    /// Analytic has its own nf_Q2 which respects the
    /// min/max nf set by lambdas
    int nf_Q2(double q2) const;
    /// Set lambda_i (for i = flavour number)
    void setLambda(unsigned int i, double lambda);

  private:

    /// Get lambdaQCD for nf
    double _lambdaQCD(int nf) const;
    /// LambdaQCD values. Stored as lambdaQCD^nf = _lambdas[nf-1]
    vector<double> _lambdas;
    /// Recalculate min/max flavors in case lambdas have changed
    void _setFlavors();
    /// Max/min number of flavors
    int _nfmax;
    int _nfmin;
  };


  /// Solve the differential equation in alphaS using an implementation of RK4
  class AlphaS_ODE : public AlphaS {
  public:
    std::string type() const { return "ode"; }
    double alphasQ2(double q2) const;
    void setMZ(double mz) { _mz = mz; }
    void setAlphaSMZ(double alphas) { _alphas_mz = alphas; }

  private:
    /// Mass of the Z-boson in GeV
    double _mz;
    /// Value of alpha_s(MZ)
    double _alphas_mz;
  };


  /// Interpolate alpha_s from tabulated points in Q2 via metadata
  class AlphaS_Ipol : public AlphaS {
  public:
    std::string type() const { return "ipol"; }
    double alphasQ2(double q2) const;
    int nf_Q2(double q2) const;

    void setQ2Values( std::vector<double> q2s ) { _q2s = q2s; }
    void setAlphaSValues( std::vector<double> as ) { _as = as; }

  private:
    double _interpolateCubic(double T, double VL, double VDL, double VH, double VDH ) const;
    double _ddq_central( size_t i ) const;
    double _ddq_forward( size_t i ) const;
    double _ddq_backward( size_t i ) const;

    std::vector<double> _q2s;
    std::vector<double> _as;
  };


}
