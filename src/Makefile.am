lib_LTLIBRARIES = libLHAPDF.la
noinst_LTLIBRARIES = libLHAPDFInfo.la libLHAPDFPaths.la libLHAPDFGlue.la

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib -avoid-version

libLHAPDFInfo_la_SOURCES = Info.cc
if USE_BUILTIN_YAMLCPP
  SUBDIRS = yamlcpp
  libLHAPDFInfo_la_CPPFLAGS = -I$(srcdir)/yamlcpp/yaml-cpp/include -DYAMLCPP_API=3 -DYAML_NAMESPACE=LHAPDF_YAML $(AM_CPPFLAGS)
  libLHAPDF_la_LDFLAGS = -L$(srcdir)/yamlcpp/yaml-cpp -llhapdf-yaml-cpp $(AM_LDFLAGS)
else
  libLHAPDFInfo_la_CPPFLAGS = $(AM_CPPFLAGS) $(YAMLCPP_CPPFLAGS) -DYAMLCPP_API=$(YAMLCPP_MAJOR_VERSION)
  libLHAPDFInfo_la_LDFLAGS = $(AM_LDFLAGS) $(YAMLCPP_LDFLAGS)
  libLHAPDFInfo_la_LIBADD = $(YAMLCPP_LIBS)
endif

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
if ENABLE_LHAGLUE
  libLHAPDF_la_LIBADD += libLHAPDFGlue.la
endif
