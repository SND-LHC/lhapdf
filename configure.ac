## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([LHAPDF],[6.0.0a0],[lhapdf@projects.hepforge.org],[LHAPDF])
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE()
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


## Library version flags (at least update the version comment with each new release version)
VERSIONFLAGS="-version-info 2:0:0" ## for 6.0.0
AC_SUBST(VERSIONFLAGS)
AC_DEFINE_UNQUOTED(LHAPDF_VERSION, "$PACKAGE_VERSION", "LHAPDF version string")


## Mac OS X compatibility
AC_CHECK_TOOL(SWVERS, sw_vers)
if test x$SWVERS != x; then
  PROD_NAME=$($SWVERS -productName | cut -f 2 -d:)
fi
AM_CONDITIONAL(WITH_OSX, [test "$PROD_NAME" == "Mac OS X"])
if test "$PROD_NAME" == "Mac OS X"; then
  MACOSX_DEPLOYMENT_TARGET=$($SWVERS -productVersion | cut -f 1,2 -d.)
  AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
  AM_CXXFLAGS="$AM_CXXFLAGS -Dunix"
  AM_CXXFLAGS="$AM_LDFLAGS -flat_namespace"
fi


## Checks and config headers
AC_LANG(C++)
AC_CONFIG_SRCDIR([src/GridPDF.cc])
AC_CONFIG_HEADERS([config/config.h include/LHAPDF/Config.h])

## Checks for programs
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL


## Enable debugging symbols.
AC_ARG_ENABLE([debug],
  [AC_HELP_STRING(--enable-debug, [build with debugging symbols (default=no)])],
  [], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_CXXFLAGS="$AM_CXXFLAGS -g"]
fi


## Enable LHAGLUE compatibility fns for Fortran
AC_ARG_ENABLE([lhaglue],
  [AC_HELP_STRING(--disable-lhaglue, [build without PDFLIB/LHAPDF5 compatibility routines])],
  [], [enable_lhaglue=yes])
if test x$enable_lhaglue == xyes; then
   AC_MSG_NOTICE(Building LHAGLUE wrappers)
else
   AC_MSG_NOTICE(Not building LHAGLUE wrappers)
fi
AM_CONDITIONAL(ENABLE_LHAGLUE, [test x$enable_lhaglue == xyes])


## Find Doxygen if possible
AC_PATH_PROG(DOXYGEN, doxygen)
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])


dnl ## Python extension
dnl AC_ARG_ENABLE(pyext, [AC_HELP_STRING(--disable-pyext,
dnl   [don't build Python module (default=build)])],
dnl   [], [enable_pyext=yes])
dnl if test x$enable_ccwrap != xyes; then
dnl    AC_MSG_WARN([Can't build Python extension since the C++ wrapper is not being built])
dnl    enable_pyext=no
dnl fi
dnl AZ_PYTHON_DEFAULT
dnl ## Basic Python checks
dnl if test x$enable_pyext == xyes; then
dnl   AZ_PYTHON_PATH
dnl   AZ_PYTHON_VERSION_ENSURE([2.3])
dnl   ## Test for Python header
dnl   if test -x "$PYTHON"; then
dnl     AC_MSG_CHECKING([for Python include path])
dnl     python_incpath=`$PYTHON -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc();"`
dnl     AC_MSG_RESULT([$python_incpath])
dnl     python_header="$python_incpath/Python.h"
dnl     if test -z "$python_incpath"; then
dnl       AC_MSG_WARN([Can't build Python extension since include directory cannot be determined from distutils])
dnl       enable_pyext=no
dnl     elif test ! -e "$python_header"; then
dnl       AC_MSG_WARN([Can't build Python extension since header file $python_header cannot be found])
dnl       enable_pyext=no
dnl     fi
dnl   fi
dnl fi
dnl ## SWIG version checks
dnl if test x$enable_pyext == xyes -a ! -e pyext/lhapdf_wrap.cc; then
dnl   AC_PATH_PROG(SWIG, swig)
dnl   if test -x "$SWIG"; then
dnl     swig_version=`$SWIG -version 2>&1 | grep Version | cut -f3 -d" "`
dnl     swig_major=`echo $swig_version | cut -f1 -d.`
dnl     swig_minor=`echo $swig_version | cut -f2 -d.`
dnl     swig_micro=`echo $swig_version | cut -f3 -d.`
dnl     let swig_version_num="1000*$swig_major+100*$swig_minor+$swig_micro"
dnl     if test "$swig_version_num" -lt 1333; then
dnl     #if test "$swig_major" -eq 1 -a "$swig_minor" -eq 3 -a "$swig_micro" -lt 33; then
dnl       AC_MSG_WARN([Can't build Python extension since SWIG version is too old to build wrapper files: please use a tarball release])
dnl       enable_pyext=no
dnl     fi
dnl   fi
dnl fi
dnl ## Finish
dnl AM_CONDITIONAL(ENABLE_PYEXT, [test x$enable_pyext == xyes])
dnl if test x$enable_pyext == xyes; then
dnl   AC_MSG_NOTICE([All Python build checks successful: 'lhapdf' Python extension will be built])
dnl fi


dnl ## Build Octave interface if possible
dnl AC_ARG_ENABLE([octave],
dnl   [AC_HELP_STRING(--disable-octave, [don't try to make Octave interface])],
dnl   [], [enable_octave=yes])
dnl if test x$enable_octave == xyes; then
dnl    AC_PATH_PROG(OCTAVE, octave)
dnl    AC_PATH_PROG(MKOCTFILE, mkoctfile)
dnl fi
dnl AM_CONDITIONAL(WITH_OCTAVE, [test x$OCTAVE != x && test x$MKOCTFILE != x])


## Boost utility library
AC_ARG_WITH([boost], AS_HELP_STRING([--with-boost], [indicate the position of the Boost library]))
if test "$with_boost"; then
  ## Check for Boost libs (only if path has been explicitly given)
  # TODO: Need lib extension independence for Macs...
  # test -e $with_boost/lib/libboost_system.so || AC_MSG_ERROR([libboost_system not found.])
  # test -e $with_boost/lib/libboost_filesystem.so || AC_MSG_ERROR([libboost_filesystem not found.])
  BOOST_CPPFLAGS=-I$with_boost/include
  BOOST_LDFLAGS=-L$with_boost/lib -R$with_boost/lib
fi
## Check for Boost headers
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS -I$prefix/include"
AC_CHECK_HEADER([boost/foreach.hpp], [], [AC_MSG_ERROR([Boost foreach not found.])])
AC_CHECK_HEADER([boost/smart_ptr.hpp], [], [AC_MSG_ERROR([Boost smart_ptr not found.])])
AC_CHECK_HEADER([boost/algorithm/string.hpp], [], [AC_MSG_ERROR([Boost algorithm/string not found.])])
AC_CHECK_HEADER([boost/lexical_cast.hpp], [], [AC_MSG_ERROR([Boost lexical_cast not found.])])
CPPFLAGS=$oldCPPFLAGS
AC_SUBST(BOOST_CPPFLAGS)
AC_SUBST(BOOST_LDFLAGS)


# yaml-cpp parser library
AC_ARG_WITH([yaml-cpp], AS_HELP_STRING([--with-yaml-cpp], [indicate the position of the yaml-cpp library]))
if test "$with_yaml_cpp"; then
  YAMLCPP_CPPFLAGS=-I$with_yaml_cpp/include
  YAMLCPP_LDFLAGS=-L$with_yaml_cpp/lib -R$with_yaml_cpp/lib
fi
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $YAMLCPP_CPPFLAGS -I$prefix/include"
AC_CHECK_HEADER([yaml-cpp/yaml.h], [], [AC_MSG_ERROR([yaml-cpp/yaml.h header not found.])])
CPPFLAGS=$oldCPPFLAGS
AC_SUBST(YAMLCPP_CPPFLAGS)
AC_SUBST(YAMLCPP_LDFLAGS)


## Set base compiler flags
AM_CPPFLAGS=["$AM_CPPFLAGS -I\$(top_srcdir)/include -I\$(prefix)/include"]
AM_CPPFLAGS=["$AM_CPPFLAGS -DLHAPDF_INSTALL_PREFIX=\\\"$prefix\\\" -DLHAPDF_DATA_PREFIX=\\\"$datadir\\\" "]
AC_SUBST(AM_CPPFLAGS)
#
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic "])
AC_CEDAR_CHECKCXXFLAG([-ansi], [AM_CXXFLAGS="$AM_CXXFLAGS -ansi "])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
AC_SUBST(AM_CXXFLAGS)


## Output
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([include/Makefile include/LHAPDF/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([doc/Makefile doc/Doxyfile])
dnl AC_CONFIG_FILES([pyext/Makefile pyext/setup.py])
dnl AC_CONFIG_FILES([octave/Makefile])
dnl AC_CONFIG_FILES([bin/Makefile bin/lhapdf-config])
AC_OUTPUT

## Message about getting PDF sets
dnl pkgdatadir=`eval "echo $datarootdir/$PACKAGE_TARNAME"`
dnl echo
dnl echo "****************************************************************"
dnl echo "IMPORTANT INFORMATION ABOUT PDF SETS"
dnl echo
dnl echo "Due to the increasing number and size of PDF data files, LHAPDF"
dnl echo "no longer bundles PDF set data in the package tarball. The sets"
dnl echo "are instead all stored online at"
dnl echo "  http://www.hepforge.org/archive/lhapdf/pdfsets/"
dnl echo "and you should install those that you wish to use into"
dnl echo "  $pkgdatadir"
dnl echo
dnl echo "The lhapdf-getdata script helps to automate the set download"
dnl echo "and installation process: before installation you will find"
dnl echo "it in the bin directory. Alternatively, download the sets by"
dnl echo "hand from the URL given above, and place them in the PDF sets"
dnl echo "directory (also given above.)"
dnl echo "****************************************************************"
dnl echo
