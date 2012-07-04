CHARMC ?= ~/workspace/charm/bin/charmc $(OPTS)
BOOST_ROOT = $(HOME)/workspace/boost_1_48_0
BOOSTINC = $(BOOST_ROOT)/include
BOOSTLIB = $(BOOST_ROOT)/lib

REVNUM  = $(shell git --git-dir=.git rev-parse HEAD)

CXX=$(CHARMC)

CXXFLAGS += -g -DAMR_REVISION=$(REVNUM) -I$(HOME)/projects/termdetection/newalg

OBJS = QuadIndex.o Advection.o Main.o

all: advection

advection: $(OBJS)
	$(CHARMC)  -module liveViz $(OPTS) $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^ -tracemode projections -balancer RefineLB

Main.decl.h Advection.decl.h: advection.ci
	$(CHARMC)  advection.ci

Advection.o: Advection.h Advection.decl.h
QuadIndex.o: 
Main.o: Main.h Advection.h Main.decl.h

png:
	g++ pngwriter.c -o pngwriter -I/home/alanger/workspace/amr/pngwriter-0.5.4/src/  -L/usr/lib -lpng -lpngwriter -lz -lfreetype



test: advection
	./charmrun ++local ./$< +p4 64 4 12

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
