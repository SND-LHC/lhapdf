lib_LTLIBRARIES = libLHAPDF.la

# PATHFINDCODE = getdatapath.cc binreloc.c

libLHAPDF_la_SOURCES = PDFGrid.cpp PDFSet.cpp
libLHAPDF_la_SOURCES += Bicubic.cpp Bilinear.cpp Nearest.cpp
libLHAPDF_la_SOURCES += AlphaS.cpp

# if ENABLE_LHAGLUE
# libLHAPDFCore_la_SOURCES += LHAGlue.cpp
# endif

if USING_WIN32
AM_CPPFLAGS += -DWIN32
else
AM_CPPFLAGS += $(BINRELOC_CFLAGS) -DENABLE_BINRELOC
endif
