lib_LTLIBRARIES = libLHAPDF.la
noinst_LTLIBRARIES = libLHAPDFInfo.la libLHAPDFPaths.la libLHAPDFGlue.la

EXTRA_DIST = yaml-cpp-0.3.0.tar.gz
if USE_BUILTIN_YAMLCPP
yaml-cpp/libyaml-cpp.a: yaml-cpp-0.3.0.tar.gz
	test -d yaml-cpp || tar xf $^; true
	grep COMPILE_FLAGS yaml-cpp/CMakeLists.txt || sed -i -e 's/LABEL_SUFFIX}"/LABEL_SUFFIX}"\n\tCOMPILE_FLAGS "-fPIC"/' yaml-cpp/CMakeLists.txt; true
	cd yaml-cpp && cmake . && $(MAKE)
all: yaml-cpp/libyaml-cpp.a
	@true
endif

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib -avoid-version

libLHAPDFInfo_la_SOURCES = Info.cc
if USE_BUILTIN_YAMLCPP
  libLHAPDFInfo_la_CPPFLAGS = -I$(srcdir)/yaml-cpp/include -DYAMLCPP_API=3 $(AM_CPPFLAGS)
  libLHAPDFInfo_la_LDFLAGS = $(AM_LDFLAGS) -Wl,-whole-archive -L$(srcdir)/yaml-cpp -lyaml-cpp -Wl,-no-whole-archive
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
