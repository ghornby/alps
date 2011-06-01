LIB_VERSION = libalps.1.2.a

CPP = g++
CINCLUDE = 


#CFLAGS = -g -DNDEBUG $(CINCLUDE)
#CFLAGS = -g -O0 -Wall -Wextra -ansi -pedantic $(CINCLUDE)
CFLAGS = -g -Wall -Wextra -ansi -pedantic $(CINCLUDE)
#CFLAGS = -O1 -Wall -Wextra -ansi -pedantic $(CINCLUDE)
#CFLAGS = -O3 -Wall $(CINCLUDE)
#CFLAGS = -O3 -Wall -finline-functions -fforce-addr -fomit-frame-pointer -fexpensive-optimizations -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 $(CINCLUDE)


## For 32 bit try either: -m32 or --target=i386


LFLAGS = -lm

#  $<  : is the name of the related file that caused the action.
#  $@  : is the name of the file to be mad.

.cpp.o:
	$(CPP) $(CFLAGS) -c $<





libalps.a:	alps.o alps_gen.o alps_sstate.o alps_utils.o alps_individ.o \
		alps_individ_real.o alps_individ_bits.o alps_layer.o \
		alps_history.o alps_random_mt.o
		ar rcv libalps.a alps.o alps_gen.o alps_sstate.o alps_utils.o \
		alps_layer.o alps_individ.o alps_individ_real.o alps_individ_bits.o \
		alps_history.o alps_random_mt.o


alps.o:		alps.h alps.cpp alps_utils.o \
		alps_individ.o alps_layer.o
		$(CPP) $(CFLAGS) -c alps.cpp

alps_gen.o:	alps.o alps_gen.h alps_gen.cpp alps_utils.o \
		alps_individ.o alps_layer.o
		$(CPP) $(CFLAGS) -c alps_gen.cpp

alps_sstate.o:	alps.o alps_sstate.h alps_sstate.cpp alps_utils.o \
		alps_individ.o alps_layer.o
		$(CPP) $(CFLAGS) -c alps_sstate.cpp

alps_layer.o:	alps_layer.h alps_layer.cpp alps_utils.o
		$(CPP) $(CFLAGS) -c alps_layer.cpp

alps_utils.o:	alps_utils.h alps_utils.cpp alps_random_mt.o
		$(CPP) $(CFLAGS) -c alps_utils.cpp

alps_individ.o:	alps_individ.h alps_individ.cpp alps_utils.o alps_history.o
		$(CPP) $(CFLAGS) -c alps_individ.cpp

alps_individ_real.o:	alps_individ_real.h alps_individ_real.cpp alps_utils.o alps_history.o
		$(CPP) $(CFLAGS) -c alps_individ_real.cpp

# Still being debugged:
alps_individ_bits.o:	alps_individ_bits.h alps_individ_bits.cpp alps_utils.o alps_history.o
		$(CPP) $(CFLAGS) -c alps_individ_bits.cpp

alps_history.o:	alps_history.h alps_history.cpp alps_utils.o alps_individ.h
		$(CPP) $(CFLAGS) -c alps_history.cpp

alps_random_mt.o:	alps_random_mt.cpp alps_random_mt.h
		$(CPP) $(CFLAGS) $(CINCLUDE) -c alps_random_mt.cpp




instlocal:	
		mkdir -p ../lib
		cp libalps.a ../lib/$(LIB_VERSION)
		rm -f ../lib/libalps.a
		ln -s $(LIB_VERSION) ../lib/libalps.a
		mkdir -p ../include/
		cp *.h ../include/

install:	
		cp libalps.a /usr/local/lib/$(LIB_VERSION)
		rm -f /usr/local/lib/libalps.a
		ln -s $(LIB_VERSION) /usr/local/lib/libalps.a
		cp *.h ../include/
		rm -r /usr/local/include/alps/
		mkdir /usr/local/include/alps/
		cp -r *.h /usr/local/include/alps/

clean:	
	rm -f *.o
	rm -f libalps.a
