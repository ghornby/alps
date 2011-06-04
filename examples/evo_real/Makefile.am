AM_CPPFLAGS = -g -Wall -I../../src -I../../include/alps
#CFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 

# We want to build evoreal but we don't want to install it.
noinst_PROGRAMS = evo_real 
evo_real_SOURCES = evo_real.cpp
evo_real_LDADD = ../../src/libalps.la
