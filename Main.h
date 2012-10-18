#ifndef MAIN_H
#define MAIN_H

#include "Main.decl.h"
class Main: public CBase_Main {

 public:
  int num_chares;
  int iterations;

  Main(CkArgMsg* m);
  //void printTreeInformation(CkVec<QuadIndex>);
  void terminate();
  void startMeshGeneration();
  //void startRunning();
  void reportCascadeStats(int *cascade_lengths, int size);
  void qdlatency(double* elems, int size);
  void remeshlatency(double* elems, int size);
  void totalWorkUnits(int total);
};

#endif
