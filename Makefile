# Charm++ compilation settings
CHARM_DIR ?= ../charm-cuda
CHARMC ?= $(CHARM_DIR)/bin/charmc -I.
CHARM_INC = -I$(CHARM_DIR)/include
CXX = $(CHARMC)
OPTS ?= -O3
DEFINE = -DTIMER
CXXFLAGS += $(DEFINE) -DAMR_REVISION=$(REVNUM) $(OPTS)

# CUDA settings
CUDATOOLKIT_HOME ?= /usr/local/cuda
NVCC ?= $(CUDATOOLKIT_HOME)/bin/nvcc
NVCC_FLAGS = -c --std=c++11 -O3
NVCC_INC = -I$(CUDATOOLKIT_HOME)/include -I$(CHARM_DIR)/src/arch/cuda/hybridAPI -I./lib/cub
CUDA_LD_LIBS = -L$(CUDATOOLKIT_HOME)/lib64 -lcudart

# Object files
SHARED_OBJS = OctIndex.o Main.o
CPU_OBJS = $(SHARED_OBJS) Advection.o
GPU_OBJS = $(SHARED_OBJS) AdvectionGPU.o AdvectionCU.o
HAPI_OBJS = $(SHARED_OBJS) AdvectionHAPI.o AdvectionHAPICU.o

# Make commands
TARGET = advection
all: $(TARGET)-cpu $(TARGET)-gpu $(TARGET)-hapi

$(TARGET)-cpu: $(CPU_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $@ $^ $(LD_LIBS) -module DistributedLB

$(TARGET)-gpu: $(GPU_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $@ $^ $(LD_LIBS) $(CUDA_LD_LIBS) -module DistributedLB

$(TARGET)-hapi: $(HAPI_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $@ $^ $(LD_LIBS) -module DistributedLB

Advection.decl.h Main.decl.h: advection.ci.stamp
advection.ci.stamp: advection.ci
	$(CHARMC) $<
	touch $@

# Compilation of object files
Advection.o: Advection.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -c $< -o $@
Main.o: Main.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -c $< -o $@
OctIndex.o: OctIndex.C OctIndex.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -c $< -o $@

AdvectionGPU.o: Advection.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -DUSE_GPU -c $< -o $@
AdvectionCU.o: Advection.cu
	$(NVCC) $(NVCC_FLAGS) $(NVCC_INC) $(CHARMINC) -o $@ $<

AdvectionHAPI.o: Advection.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -DUSE_GPU -DUSE_HAPI -c $< -o $@
AdvectionHAPICU.o: Advection.cu
	$(NVCC) $(NVCC_FLAGS) -DUSE_HAPI $(NVCC_INC) $(CHARMINC) -o $@ $<

# Tests are currently set for SMP
test: $(TARGET)-cpu
	./charmrun +p8 ./$(TARGET)-cpu 3 32 30 9 +balancer DistributedLB ++ppn 8

test-gpu: $(TARGET)-gpu
	./charmrun +p8 ./$(TARGET)-gpu 3 32 30 9 +balancer DistributedLB ++ppn 8

test-hapi: $(TARGET)-hapi
	./charmrun +p8 ./$(TARGET)-hapi 3 32 30 9 +balancer DistributedLB ++ppn 8

clean:
	rm -f *.decl.h *.def.h conv-host *.o $(TARGET)-cpu $(TARGET)-gpu $(TARGET)-hapi charmrun advection.ci.stamp
