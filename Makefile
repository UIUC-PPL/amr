
CHARMHOME ?= ~/programming/ppl/charm/net-darwin-x86_64
CHARMC ?= $(CHARMHOME)/bin/charmc -I.
CXX=$(CHARMC)

OPTS ?= -O0 -g
CXXFLAGS += -DAMR_REVISION=$(REVNUM) $(OPTS)

OBJS = OctIndex.o Advection.o Main.o 

all: advection

advection: $(OBJS)
	$(CHARMC) $(CXXFLAGS) $(LDFLAGS) -language charm++ -o $@ $^

Advection.decl.h Main.decl.h: advection.ci.stamp
advection.ci.stamp: advection.ci
	$(CHARMC) $<
	touch $@

Advection.o: Advection.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
Main.o: Main.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
OctIndex.o: OctIndex.C OctIndex.h Advection.decl.h

debug: all
	~/programming/ppl/ccs_tools/bin/charmdebug ./advection 3 32 30 9 +p4

test: all
	./charmrun +p4 ++local ./advection 5 32 30 9

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun advection.ci.stamp log/*log

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
