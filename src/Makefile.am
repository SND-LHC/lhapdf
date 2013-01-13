lib_LTLIBRARIES = libLHAPDF.la

libLHAPDF_la_SOURCES = \
  PDFGrid.cc Info.cc Paths.cc \
  BicubicInterpolator.cc BilinearInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc Factories.cc

if ENABLE_LHAGLUE
  libLHAPDF_la_SOURCES += LHAGlue.cc
endif

libLHAPDF_la_LIBADD = -lboost_filesystem -lboost_system -lyaml-cpp
