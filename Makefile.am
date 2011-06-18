AUTOMAKE_OPTIONS = foreign
SUBDIRS = src examples support include/alps
ACLOCAL_AMFLAGS = -I support/m4

# include support/aminclude.am

# if DX_COND_man

# # You'd probably want to post-process man pages and installed the patched
# # versions.
# #dist_man1_MANS = @DX_DOCDIR@/man/man1/doxample.c.1

# $(man3_MANS): doxygen-doc

# endif

# doc: doxygen-doc

# html-local: doxygen-doc

MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = README.txt LICENSE.txt $(DX_CONFIG) doc/api/html

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = support/alps.pc
