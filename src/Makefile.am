lib_LTLIBRARIES = libLHAPDF.la

# PATHFINDCODE = getdatapath.cc binreloc.c

libLHAPDF_la_SOURCES = PDFGrid.cc PDFSet.cc
libLHAPDF_la_SOURCES += Bicubic.cc Bilinear.cc Nearest.cc
libLHAPDF_la_SOURCES += AlphaS.cc

# if ENABLE_LHAGLUE
# libLHAPDFCore_la_SOURCES += LHAGlue.cc
# endif

if USING_WIN32
AM_CPPFLAGS += -DWIN32
else
AM_CPPFLAGS += $(BINRELOC_CFLAGS) -DENABLE_BINRELOC
endif
