lib_LTLIBRARIES = libLHAPDF.la

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS)
AM_CPPFLAGS += -DLHAPDF_INSTALL_PREFIX=\"$(prefix)\" -DLHAPDF_DATA_PREFIX=\"$(datadir)\" -DYAMLCPP_API=$(YAMLCPP_MAJOR_VERSION)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS)

libLHAPDF_la_SOURCES = \
  GridPDF.cc Paths.cc \
  Info.cc PDFInfo.cc \
  Interpolator.cc BilinearInterpolator.cc BicubicInterpolator.cc \
  LogBilinearInterpolator.cc LogBicubicInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc AlphaS_Analytic.cc AlphaS_ODE.cc AlphaS_Ipol.cc \
  Factories.cc

if ENABLE_LHAGLUE
  libLHAPDF_la_SOURCES += LHAGlue.cc
endif

libLHAPDF_la_LIBADD = $(BOOST_LIBS) $(YAMLCPP_LIBS)

# libLHAPDF_la_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS) -DYAMLCPP_API=$(YAMLCPP_MAJOR_VERSION)
# libLHAPDF_la_LDFLAGS = $(AM_LDFLAGS) -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS) -export-dynamic
