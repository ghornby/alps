#AM_CPPFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 
AM_CPPFLAGS =-Wall -Wextra -ansi -pedantic -I../include/alps

lib_LTLIBRARIES = libalps.la

# The alps_individ_bits class not included; it breaks the build.
libalps_la_SOURCES = alps.cpp alps_gen.cpp alps_history.cpp alps_individ.cpp alps_individ_real.cpp alps_layer.cpp alps_random_mt.cpp alps_sstate.cpp alps_utils.cpp

# There are two library versioning systems.  One aligns with the
# release number.  The other is the version-info which has a trickier
# specification.  Below is the relevant gory details.  We're going to
# just use the release number.  -Shane

# Updating library version information
# ====================================
#
#    If you want to use libtool's versioning system, then you must specify
# the version information to libtool using the `-version-info' flag
# during link mode (*note Link mode::).
#
#    This flag accepts an argument of the form
# `CURRENT[:REVISION[:AGE]]'.  So, passing `-version-info 3:12:1' sets
# CURRENT to 3, REVISION to 12, and AGE to 1.
#
#    If either REVISION or AGE are omitted, they default to 0.  Also note
# that AGE must be less than or equal to the CURRENT interface number.
#
#    Here are a set of rules to help you update your library version
# information:
#
#   1. Start with version information of `0:0:0' for each libtool library.
#
#   2. Update the version information only immediately before a public
#      release of your software.  More frequent updates are unnecessary,
#      and only guarantee that the current interface number gets larger
#      faster.
#
#   3. If the library source code has changed at all since the last
#      update, then increment REVISION (`C:R:A' becomes `C:r+1:A').
#
#   4. If any interfaces have been added, removed, or changed since the
#      last update, increment CURRENT, and set REVISION to 0.
#
#   5. If any interfaces have been added since the last public release,
#      then increment AGE.
#
#   6. If any interfaces have been removed since the last public release,
#      then set AGE to 0.
#
#    *_Never_* try to set the interface numbers so that they correspond
# to the release number of your package.  This is an abuse that only
# fosters misunderstanding of the purpose of library versions.  Instead,
# use the `-release' flag (*note Release numbers::), but be warned that
# every release of your package will not be binary compatible with any
# other release.
#libalps_la_LDFLAGS = -version-info 0:0:0
libalps_la_LDFLAGS = -release $(PACKAGE_VERSION)


