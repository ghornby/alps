AM_CPPFLAGS = -Wall -Wextra -ansi -pedantic -I../../src
#CFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 

# We want to build evoreal but we don't want to install it.
noinst_PROGRAMS =  evo_real_barebones

evo_real_barebones_SOURCES = evo_real_barebones.cpp
evo_real_barebones_LDADD = ../../src/libalps.la
