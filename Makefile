
CHARMHOME ?= $(HOME)/git/charm/net-linux-x86_64-syncft
CHARMC ?= $(CHARMHOME)/bin/charmc -module liveViz
CXX=$(CHARMC)

OPTS ?= -O3
CXXFLAGS += -g -DLOGGER -DAMR_REVISION=$(REVNUM) $(OPTS)

OBJS = QuadIndex.o Advection.o Main.o

all: advection

advection: $(OBJS)
	$(CHARMC) $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^ -balancer RotateLB -lboost_filesystem -lboost_system

Advection.decl.h Main.decl.h: advection.ci.stamp
advection.ci.stamp: advection.ci
	$(CHARMC) $<
	touch $@

Advection.o: Advection.C Advection.h QuadIndex.h Main.decl.h Advection.decl.h
Main.o: Main.C Advection.h QuadIndex.h Main.decl.h Advection.decl.h
QuadIndex.o: QuadIndex.C QuadIndex.h Advection.decl.h

test: all
	./charmrun +p4 ./advection 5 16 20

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun advection.ci.stamp

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
