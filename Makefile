CHARMC=~/charm/bin/charmc -module liveViz $(OPTS)
CXX=$(CHARMC)
OBJS = QuadIndex.o Advection.o

all: advection

advection: $(OBJS)
	$(CHARMC) $(OPTS) -language charm++ -o advection Main.C  $(OBJS)

advection.decl.h: advection.ci
	$(CHARMC)  advection.ci

Advection.o: advection.decl.h
	$(CHARMC) $(OPTS) -o Advection.o Advection.C

QuadIndex.o: 
	$(CHARMC) $(OPTS) -o QuadIndex.o QuadIndex.C

test: all
	./charmrun advetion +p4 10

clean:
	rm -f *.decl.h *.def.h conv-host *.o advection charmrun

bgtest: all
	./charmrun advection +p4 10 +x2 +y2 +z2
