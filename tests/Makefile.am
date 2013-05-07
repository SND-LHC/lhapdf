FLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) $(YAMLCPP_CPPFLAGS) $(AM_LDFLAGS) $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS) -L$(top_builddir)/src -lLHAPDF $(BOOST_LIBS) $(YAMLCPP_LIBS)

.PHONY = all testalphas testgrid testindex testinfo testpaths

all: testalphas testgrid testindex testinfo testpaths
	@true

testalphas: testalphas.cc
	$(CXX) testalphas.cc $(FLAGS) -o testalphas

testgrid: testgrid.cc
	$(CXX) testgrid.cc $(FLAGS) -o testgrid

testindex: testindex.cc
	$(CXX) testindex.cc $(FLAGS) -o testindex

testinfo: testinfo.cc
	$(CXX) testinfo.cc $(FLAGS) -o testinfo

testpaths: testpaths.cc
	$(CXX) testpaths.cc $(FLAGS) -o testpaths
