noinst_PROGRAMS = testpdf testpdfset analyticpdf compatibility

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib $(BOOST_LDFLAGS)
LIBS = -lLHAPDF $(BOOST_LIBS)

testpdf_SOURCES = testpdf.cc
testpdfset_SOURCES = testpdfset.cc
analyticpdf_SOURCES = analytic.cc
compatibility_SOURCES = compatibility.cc
