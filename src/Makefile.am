SUBDIRS = yamlcpp

lib_LTLIBRARIES = libLHAPDF.la
noinst_LTLIBRARIES = libLHAPDFInfo.la libLHAPDFPaths.la libLHAPDFGlue.la

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib -avoid-version

libLHAPDFInfo_la_SOURCES = Info.cc
libLHAPDFInfo_la_CPPFLAGS = -I$(srcdir)/yamlcpp/yaml-cpp/include -DYAMLCPP_API=3 -DYAML_NAMESPACE=LHAPDF_YAML $(AM_CPPFLAGS)

libLHAPDFPaths_la_SOURCES = Paths.cc
libLHAPDFPaths_la_CPPFLAGS = $(AM_CPPFLAGS) -DLHAPDF_INSTALL_PREFIX=\"$(prefix)\" -DLHAPDF_DATA_PREFIX=\"$(datadir)\"

libLHAPDFGlue_la_SOURCES = LHAGlue.cc
if ENABLE_LHAGLUE_CXX
libLHAPDFGlue_la_CPPFLAGS = $(AM_CPPFLAGS) -DENABLE_LHAGLUE_CXX
endif

libLHAPDF_la_SOURCES = \
  PDF.cc GridPDF.cc PDFInfo.cc \
  Interpolator.cc BilinearInterpolator.cc BicubicInterpolator.cc \
  LogBilinearInterpolator.cc LogBicubicInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc AlphaS_Analytic.cc AlphaS_ODE.cc AlphaS_Ipol.cc \
  Factories.cc

libLHAPDF_la_LIBADD = $(BOOST_LIBS) libLHAPDFInfo.la libLHAPDFPaths.la
libLHAPDF_la_LIBADD += $(srcdir)/yamlcpp/yaml-cpp/liblhapdf-yaml-cpp.a
if ENABLE_LHAGLUE
  libLHAPDF_la_LIBADD += libLHAPDFGlue.la
endif
libLHAPDF_la_LIBTOOLFLAGS = --silent
