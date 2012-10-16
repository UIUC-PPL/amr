CHARMC=~/charm/bin/charmc -module liveViz -DLOGGER
BOOST_ROOT = $(HOME)/boost_1_46_1
BOOSTINC = $(BOOST_ROOT)/include
BOOSTLIB = $(BOOST_ROOT)/lib

CXX=$(CHARMC)

OPTS ?= -O3
CXXFLAGS += -g -DAMR_REVISION=$(REVNUM) $(OPTS)

OBJS = QuadIndex.o Advection.o Main.o

all: advection

advection: $(OBJS)
	$(CHARMC)  -module liveViz $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^ -tracemode projections -balancer RotateLB -lboost_filesystem

Advection.decl.h Main.decl.h: advection.ci.stamp
advection.ci.stamp: advection.ci
	$(CHARMC) $<
	touch $@

Advection.o: Advection.C Advection.h Main.decl.h Advection.decl.h
Main.o: Main.C Advection.h Main.decl.h Advection.decl.h
QuadIndex.o: QuadIndex.C QuadIndex.h

test: all
	./charmrun advetion +p4 10

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun advection.ci.stamp

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
