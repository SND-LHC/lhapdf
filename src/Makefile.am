lib_LTLIBRARIES = libLHAPDF.la

# PATHFINDCODE = getdatapath.cc binreloc.c

libLHAPDF_la_SOURCES = \
  PDFGrid.cc \
  BicubicInterpolator.cc BilinearInterpolator.cc \
  ErrExtrapolator.cc NearestPointExtrapolator.cc
# AlphaS.cc
#PDFSet.cc

# if ENABLE_LHAGLUE
# libLHAPDFCore_la_SOURCES += LHAGlue.cc
# endif

if USING_WIN32
AM_CPPFLAGS += -DWIN32
else
AM_CPPFLAGS += $(BINRELOC_CFLAGS) -DENABLE_BINRELOC
endif
