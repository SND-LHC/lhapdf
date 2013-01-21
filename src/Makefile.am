lib_LTLIBRARIES = libLHAPDF.la

libLHAPDF_la_SOURCES = \
  GridPDF.cc Info.cc Paths.cc \
  BicubicInterpolator.cc BilinearInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc AlphaS_Analytic.cc AlphaS_ODE.cc \
  Factories.cc

if ENABLE_LHAGLUE
  libLHAPDF_la_SOURCES += LHAGlue.cc
endif

BOOSTDIR=/usr # TODO: determine this and YAMLDIR in configure.ac
libLHAPDF_la_LDFLAGS = -L$(prefix)/lib -L$(BOOSTDIR)/lib -export-dynamic
libLHAPDF_la_LIBADD = -lboost_filesystem -lboost_system -lyaml-cpp
# TODO: How to make the Boost libs auto-chain to libs linking against LHAPDF?
