# Charm++ compilation settings
CHARM_DIR ?= ../charm
CHARMC ?= $(CHARM_DIR)/bin/charmc -I. $(CHARMC_INC)
CHARMC_INC = -I$(CUDATOOLKIT_HOME)/include -I$(CHARM_DIR)/include
CXX = $(CHARMC)
OPTS ?= --std=c++11 -O3
CXXFLAGS += $(OPTS)

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
all: cpu gpu hapi

cpu: $(CPU_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $(TARGET)-$@ $^ $(LD_LIBS) -module CommonLBs

gpu: $(GPU_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $(TARGET)-$@ $^ $(LD_LIBS) $(CUDA_LD_LIBS) -module CommonLBs

hapi: $(HAPI_OBJS)
	$(CHARMC) $(CXXFLAGS) -language charm++ -o $(TARGET)-$@ $^ $(LD_LIBS) $(CUDA_LD_LIBS) -module CommonLBs

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
	$(NVCC) $(NVCC_FLAGS) $(NVCC_INC) -o $@ $<

AdvectionHAPI.o: Advection.C Advection.h OctIndex.h Main.decl.h Advection.decl.h
	$(CHARMC) $(CXXFLAGS) -DUSE_GPU -DUSE_HAPI -c $< -o $@
AdvectionHAPICU.o: Advection.cu
	$(NVCC) $(NVCC_FLAGS) -DUSE_HAPI $(NVCC_INC) -o $@ $<

# Tests are currently set for SMP
test: cpu
	./charmrun +p8 ./$(TARGET)-$< -a 128 -b 32 -d 3 -i 30 -l 9 +balancer DistributedLB ++ppn 8 ++local

test-gpu: gpu
	./charmrun +p8 ./$(TARGET)-$< -a 128 -b 32 -d 3 -i 30 -l 9 +balancer DistributedLB ++ppn 8 ++local

test-hapi: hapi
	./charmrun +p8 ./$(TARGET)-$< -a 128 -b 32 -d 3 -i 30 -l 9 +balancer DistributedLB ++ppn 8 ++local

clean:
	rm -f *.decl.h *.def.h conv-host *.o $(TARGET)-cpu $(TARGET)-gpu $(TARGET)-hapi charmrun advection.ci.stamp
