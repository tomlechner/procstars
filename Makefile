#################################################
##############                    ###############
#############  Tycho Generator     ##############
##############                    ###############
#################################################


LAXDIR=/home/tom/other/p/sourceforge/laxkit/lax

BINDIR=$(PREFIX)/bin


LD=g++
CPPFLAGS= -Wall -g -I. -I.. -I$(LAXDIR)/.. -I$(LAXIDIR) -I/usr/include/freetype2 -I/usr/include/GraphicsMagick/
LDFLAGS= -L/usr/X11R6/lib -lX11 -lftgl -lm -lpng -lcrypto -lGraphicsMagick++ `imlib2-config --libs` -L$(LAXDIR) -L$(LAXIDIR) -lXft -lXi






#tychomaker: lax
process-tycho: process-tycho.cc
	g++ $(pobjs) process-tycho.cc $(LDFLAGS) $(CPPFLAGS) -o $@
	#g++ $(pobjs) process-tycho.o -llaxinterfaces -llaxkit $(TUIOLIBS) $(LDFLAGS) -lGL -lGLU $(CPPFLAGS) -o $@

test: test.cc
	g++ $(pobjs) test.cc $(LDFLAGS) $(CPPFLAGS) -o $@

lax:
	cd $(LAXDIR) && $(MAKE)
	cd $(LAXDIR)/interfaces && $(MAKE)

laxinterface:
	cd $(LAXDIR)/interfaces && $(MAKE)


#depends:
#	makedepend -fmakedepend -I$(LAXDIR)/.. -Y *.cc
#
#include makedepend




.PHONY: clean lax laxinterface
clean:
	rm -f tycho *.o

