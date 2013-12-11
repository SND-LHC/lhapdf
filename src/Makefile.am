lib_LTLIBRARIES = libLHAPDF.la
noinst_LTLIBRARIES = libLHAPDFInfo.la libLHAPDFPaths.la libLHAPDFGlue.la

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) -avoid-version

libLHAPDFInfo_la_SOURCES = Info.cc
libLHAPDFInfo_la_CPPFLAGS = $(AM_CPPFLAGS) $(YAMLCPP_CPPFLAGS) -DYAMLCPP_API=$(YAMLCPP_MAJOR_VERSION)
libLHAPDFInfo_la_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS) -avoid-version
libLHAPDFInfo_la_LIBADD =  $(YAMLCPP_LIBS)

libLHAPDFPaths_la_SOURCES = Paths.cc
libLHAPDFPaths_la_CPPFLAGS = $(AM_CPPFLAGS) -DLHAPDF_INSTALL_PREFIX=\"$(prefix)\" -DLHAPDF_DATA_PREFIX=\"$(datadir)\"
libLHAPDFPaths_la_LIBADD = $(BOOST_LIBS)

libLHAPDFGlue_la_SOURCES = LHAGlue.cc

libLHAPDF_la_SOURCES = \
  PDF.cc GridPDF.cc PDFInfo.cc \
  Interpolator.cc BilinearInterpolator.cc BicubicInterpolator.cc \
  LogBilinearInterpolator.cc LogBicubicInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc \
  AlphaS.cc AlphaS_Analytic.cc AlphaS_ODE.cc AlphaS_Ipol.cc \
  Factories.cc

libLHAPDF_la_LIBADD = $(BOOST_LIBS) libLHAPDFInfo.la
if ENABLE_LHAGLUE
  libLHAPDF_la_LIBADD += libLHAPDFGlue.la
endif
