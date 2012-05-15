CHARMC=~/projects/termdetection/charm/bin/charmc -module liveViz $(OPTS)
BOOST_ROOT = $(HOME)/boost
BOOSTINC = $(BOOST_ROOT)/include
BOOSTLIB = $(BOOST_ROOT)/lib

REVNUM  = $(shell git --git-dir=.git rev-parse HEAD)

CXX=$(CHARMC)

CXXFLAGS += -I$(BOOSTINC) -g -DAMR_REVISION=$(REVNUM) -I$(HOME)/projects/termdetection/newalg
LDFLAGS += -L$(BOOSTLIB)

OBJS = QuadIndex.o Advection.o Main.o

all: advection

advection: $(OBJS)
	$(CHARMC) $(OPTS) $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^ -lboost_filesystem -tracemode projections

Main.decl.h Advection.decl.h: advection.ci
	$(CHARMC)  advection.ci

Advection.o: Advection.h Advection.decl.h
QuadIndex.o: 
Main.o: Main.h Advection.h Main.decl.h

test: advection
	./charmrun ++local ./$< +p4 64 4 12

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
