lib_LTLIBRARIES = libLHAPDF.la

libLHAPDF_la_SOURCES = \
  GridPDF.cc Info.cc Paths.cc \
  Interpolator.cc BilinearInterpolator.cc BicubicInterpolator.cc \
  LogBilinearInterpolator.cc LogBicubicInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc AlphaS_Analytic.cc AlphaS_ODE.cc AlphaS_Ipol.cc \
  Factories.cc

if ENABLE_LHAGLUE
  libLHAPDF_la_SOURCES += LHAGlue.cc
endif

libLHAPDF_la_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS)
libLHAPDF_la_LDFLAGS = $(AM_LDFLAGS) -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS) -export-dynamic
libLHAPDF_la_LIBADD = $(BOOST_LIBS) -lyaml-cpp
