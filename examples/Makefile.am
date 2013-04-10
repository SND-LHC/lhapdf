bin_PROGRAMS = analyticpdf compatibility

AM_CPPFLAGS += $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS)
AM_LDFLAGS = -L$(top_srcdir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS)
LIBS = -lLHAPDF $(BOOST_LIBS) $(YAMLCPP_LIBS)

analyticpdf_SOURCES = analytic.cc
compatibility_SOURCES = compatibility.cc
