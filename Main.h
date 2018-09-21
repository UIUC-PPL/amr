#ifndef MAIN_H_
#define MAIN_H_

#include "Main.decl.h"
class Main : public CBase_Main {
 public:
  int num_chares;
  int iterations;
  int nresponses;

  Main(CkArgMsg* m);
  Main(CkMigrateMessage* m) : CBase_Main(m) {}
  void terminate();
  void startMeshGeneration();
  void totalWorkUnits(int total);
};

#endif // MAIN_H_
