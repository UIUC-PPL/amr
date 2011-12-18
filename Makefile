CHARMC=~/charm/bin/charmc -module liveViz $(OPTS)
BOOST_ROOT = $(HOME)/boost
BOOSTINC = $(BOOST_ROOT)/include
BOOSTLIB = $(BOOST_ROOT)/lib

REVNUM  = $(shell git --git-dir=.git rev-parse HEAD)

CXX=$(CHARMC)

CPPFLAGS += -I$(BOOSTINC) -DAMR_REVISION=$(REVNUM)
LDFLAGS += -L$(BOOSTLIB)

OBJS = QuadIndex.o Advection.o

all: advection

advection: $(OBJS)
	$(CHARMC) $(OPTS) $(CPPFLAGS) $(LDFLAGS) -language charm++ -o advection Main.C  $(OBJS) -lboost_filesystem -tracemode projections

advection.decl.h: advection.ci
	$(CHARMC)  advection.ci

Advection.o: advection.decl.h
	$(CHARMC) $(OPTS) $(CPPFLAGS) $(LDFLAGS) -o Advection.o Advection.C

QuadIndex.o: 
	$(CHARMC) $(OPTS) $(CPPFLAGS) $(LDFLAGS) -o QuadIndex.o QuadIndex.C

test: all
	./charmrun advetion +p4 10

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
