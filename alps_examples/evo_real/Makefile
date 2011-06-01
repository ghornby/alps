CINCLUDE = -I../include

CFLAGS = -g -Wall $(CINCLUDE)
#CFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 $(CINCLUDE)

CC = gcc



LIBSUFFIX = .a

LINKALPS = -L../lib -lalps

LINK = $(LINKALPS) -lm -lstdc++

TARGETS = \
	evoreal


all:	$(TARGETS)


clean:
	$(RM) *.o
	$(RM) $(TARGETS)


evoreal.o:	evoreal.cpp
		$(CC) $(CFLAGS) -c evoreal.cpp

evoreal:	evoreal.o
		$(CC) $(CFLAGS) evoreal.o $(LINK) -o evoreal

