FLAGS := $(AM_CPPFLAGS) -I$(prefix)/include $(AM_LDFLAGS) -L$(prefix)/lib -L$(top_srcdir)/src -lLHAPDF -lboost_filesystem -lboost_system

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
