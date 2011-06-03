#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT([alps], [1.1.2], [bug-report@address])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([support/m4])

# Document generation
dnl DX_HTML_FEATURE(ON)
dnl DX_CHM_FEATURE(OFF)
dnl DX_CHI_FEATURE(OFF)
dnl DX_MAN_FEATURE(OFF)
dnl DX_RTF_FEATURE(OFF)
dnl DX_XML_FEATURE(OFF)
dnl DX_PDF_FEATURE(OFF)
dnl DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(alps, doc/doxygen.cfg, doc/api)


AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m], [fabs])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor gettimeofday pow sqrt])

# Output.
AM_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([Makefile 
                 examples/Makefile
                 examples/evo_real/Makefile
                 examples/evo_real_barebones/Makefile
                 src/Makefile
                 support/Makefile
                 support/alps-config
                 support/alps.pc])
AC_OUTPUT
