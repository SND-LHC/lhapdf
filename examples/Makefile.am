bin_PROGRAMS = testpdf testpdfset analyticpdf compatibility

AM_CPPFLAGS += $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS)
AM_LDFLAGS = -L$(top_srcdir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS)
LIBS = -lLHAPDF $(BOOST_LIBS) $(YAMLCPP_LIBS)

testpdf_SOURCES = testpdf.cc
testpdfset_SOURCES = testpdfset.cc
analyticpdf_SOURCES = analytic.cc
compatibility_SOURCES = compatibility.cc

# TODO: need to remove this, but it was fun
testpdfset_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11
