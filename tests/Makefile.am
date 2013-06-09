noinst_PROGRAMS = testalphas testgrid testindex testinfo testpaths testperf

AM_CPPFLAGS += -I$(top_srcdir)/include -I$(prefix)/include $(BOOST_CPPFLAGS)
AM_LDFLAGS += -L$(top_builddir)/src -L$(prefix)/lib $(BOOST_LDFLAGS)
LIBS = -lLHAPDF $(BOOST_LIBS)

testalphas_SOURCES = testalphas.cc
testgrid_SOURCES = testgrid.cc
testindex_SOURCES = testindex.cc
testinfo_SOURCES = testinfo.cc
testpaths_SOURCES = testpaths.cc
testperf_SOURCES = testperf.cc
