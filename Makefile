
CHARMHOME ?= $(HOME)/charm
CHARMC ?= $(CHARMHOME)/bin/charmc
CXX=$(CHARMC)

OPTS ?= -O3
CXXFLAGS += -g -DAMR_REVISION=$(REVNUM) $(OPTS)

OBJS = QuadIndex.o Advection.o Main.o

all: advection

advection: $(OBJS)
	$(CHARMC) $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^ -tracemode projections -balancer RotateLB 

Advection.decl.h Main.decl.h: advection.ci.stamp
advection.ci.stamp: advection.ci
	$(CHARMC) $<
	touch $@

Advection.o: Advection.C Advection.h QuadIndex.h Main.decl.h Advection.decl.h
Main.o: Main.C Advection.h QuadIndex.h Main.decl.h Advection.decl.h
QuadIndex.o: QuadIndex.C QuadIndex.h Advection.decl.h

test: all
	./charmrun +p8 ++local ./advection 5 16 4

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun advection.ci.stamp

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
