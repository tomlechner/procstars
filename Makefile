#################################################
##############                    ###############
#############      Procstars       ##############
##############                    ###############
#################################################


#--for both command line only, and gui versions:
GRAPHICSMAGICKDIR=/usr/include/GraphicsMagick/


#---for gui version:
FREETYPEDIR=/usr/include/freetype2
LAXDIR=./laxkit/lax




#--------- hopefully you don't have to mess with anything below ----------


LD=g++

LAXCPPFLAGS= -Wall -g -I. -I.. -I$(LAXDIR)/.. -I$(FREETYPEDIR) 
LAXLDFLAGS= -L/usr/X11R6/lib -lX11 -lftgl -lm -lpng -lcrypto -lcairo `imlib2-config --libs`  -lXft -lXi -L$(LAXDIR) 


CPPFLAGS= $(LAXCPPFLAGS) -Wall -g -I.  -I$(GRAPHICSMAGICKDIR)
LDFLAGS= -L/usr/X11R6/lib -lm  -lGraphicsMagick++ 





procstars: procstars.cc
	g++ procstars.cc $(LDFLAGS) $(CPPFLAGS) -o $@

galaxymaker: galaxymaker.cc
	g++ $@.cc $(LDFLAGS) $(CPPFLAGS) -o $@

procstars-gui: lax catalogs.o procstars-gui.o
	g++  $(LDFLAGS) $(LAXLDFLAGS) $@.o catalogs.o -llaxkit $(CPPFLAGS)  -o $@

test: test.cc
	g++  test.cc $(LDFLAGS) $(CPPFLAGS) -o $@

lax:
	cd $(LAXDIR) && $(MAKE)
	cd $(LAXDIR)/interfaces && $(MAKE)

laxinterface:
	cd $(LAXDIR)/interfaces && $(MAKE)


depends:
	makedepend -fmakedepend -I$(LAXDIR)/.. -Y *.cc

include makedepend




.PHONY: clean lax laxinterface depends
clean:
	rm -f procstars-gui procstars *.o

