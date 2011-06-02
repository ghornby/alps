AM_CPPFLAGS = -g -Wall -I../../alps_src
#CFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 

# We want to build evoreal but we don't want to install it.
noinst_PROGRAMS = evoreal 
evoreal_SOURCES = evoreal.cpp
evoreal_LDADD = ../../alps_src/libalps.la
