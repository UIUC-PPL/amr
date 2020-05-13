#ifndef ADVECTION_H_
#define ADVECTION_H_

#include "Advection.decl.h"

#ifdef USE_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#endif

class Neighbor {
  bool refined;
  bool dataSent;
  Decision decision;
  Decision childDecisions[NUM_CHILDREN];
  int dir;

public:
  Neighbor() : refined(false), dataSent(false), decision(INV), dir(-1) {}
  Neighbor(int dir) : refined(false), dataSent(false), decision(INV), dir(dir) {}

  int getDir() { return dir; }

  bool isRefined() {
    return refined;
  }

  bool isDataSent() {
    return dataSent;
  }

  void setDataSent(bool value) {
    dataSent = value;
  }

  void setRefined(bool value) {
    refined = value;
  }

  Decision getDecision(int child = -1) {
    assert((child == -1) == !refined);
    Decision &D = (child == -1) ? decision : childDecisions[child];
    return D;
  }

  Decision setDecision(Decision d, int child = -1) {
    assert((child == -1) == !refined);
    Decision &D = (child == -1) ? decision : childDecisions[child];
    return D = std::max(D, d);
  }

  void pup(PUP::er &p) {
    p|refined;
    p|dataSent;
    p|decision;
    p|dir;
    PUParray(p, childDecisions, NUM_CHILDREN);
  }

  void resetDecision() {
    decision = COARSEN;
    for (int i = 0; i < NUM_CHILDREN; ++i)
      childDecisions[i] = COARSEN;
  }
};

class MeshBlock : public CBase_MeshBlock {
  MeshBlock_SDAG_CODE
public:
  MeshManager* mesh_manager_local;

  // Tree information
  bool isRefined;
  int depth;
  bool child_isRefined[NUM_CHILDREN];
  bool isGrandParent();

  std::set<OctIndex> ghostReceived;

  // Phase 1 data structures
  Decision decision;
  Decision child_decision[NUM_CHILDREN];

  bool parentHasAlreadyMadeDecision;//to be used by a parent
  bool hasReceivedParentDecision;
  char fname[100];
  OctIndex  parent;
  std::map<OctIndex, Neighbor> neighbors;
  std::map<OctIndex, Decision> uncleDecisions;
  int xc, yc, zc;

  //data
  float imsg;

  float* u;
  float* u2;
  float* u3;
  float *x;
  float *y;
  float *z;

#ifdef USE_GPU
  float* d_u;
  float* d_u2;
  float* d_u3;
  float* d_error;
  float* h_error;

  cudaStream_t computeStream;
  cudaStream_t decisionStream;
#endif

  float *left_surface;
  float *right_surface;
  float *top_surface;
  float *bottom_surface;
  float *forward_surface;
  float *backward_surface;

  int iterations;
  int meshGenIterations;

  float up_x, up_y, up_z;
  float un_x, un_y, un_z;
  float myt, mydt;
  float dx, dy, dz, nx, ny, nz;
  float x_min, x_max, y_min, y_max, z_min, z_max;
  
  float itBeginTime;
  float lastBusyTime, lastIdleTime;
  float remeshStartTime, remeshEndTime;
  void mem_allocate(float* &p, int size);
  void mem_allocate_all();
  void mem_deallocate_all();
  OctIndex getRefinedNeighbor(int NBR);
  int getSourceDirection(int NBR);
  float* getGhostBuffer(int dir);
  int getGhostCount(int dir);
  float lastIdleTimeQD;
  std::ofstream logfile;

  bool finishedPhase1;
  int nChildDataRecvd;
  bool phase1Over;

  double compute_start_time;
  double decision_start_time;

  ~MeshBlock();

  /* Constructors */
  MeshBlock() { usesAtSync = true; usesAutoMeasure = false; }
  MeshBlock(CkMigrateMessage* m) : CBase_MeshBlock(m){ 
    usesAutoMeasure = false;
    usesAtSync = true;
    VB(logfile.open(string("log/"+thisIndex.getIndexString()+"log").c_str(), std::ofstream::app););
    VB(logfile << "migrated to a new processor" << std::endl;)
  }
  MeshBlock(float, float, float, float, float, float);
  MeshBlock(float dx, float dy, float dz,
            float myt, float mydt,
            float x_min, float y_min, float z_min,
            int meshGenIterations_, int iterations_,
            std::vector<float> refined_u,
            std::map<OctIndex, Neighbor> neighbors);

  void initializeRestofTheData(); // common function initialization for
  
  void pup(PUP::er &p);

  /* initial mesh generation*/
  void applyInitialCondition();
  void process(int, int, int, int, float*);
  void compute();
  void computeDone(); // GPU - HAPI

  /*Phase1 Entry Methods*/
  void makeGranularityDecisionAndCommunicate();
  Decision getGranularityDecision();
  void gotErrorFromGPU(); // GPU - HAPI

  void resetMeshRestructureData();
  void prepareData4Exchange();
  void processPhase1Msg(int, int, Decision, int);

  void updateDecisionState(int cascade_length, Decision newDecision);

  void processChildDecision( int, Decision, int);

  void processParentDecision(int cascade_length);

  /*Phase2 entry methods*/
  void setNbrStateUponCoarsening(int, int, std::map<OctIndex, Neighbor>&, std::map<OctIndex, Decision> &);
  void sendReadyData();
  // Returns whether a message was sent
  void sendGhost(int);
  void doPhase2();
  void updateMeshState();
  
  void recvChildData(int, float, float, int, int, std::vector<float>, std::map<OctIndex, Neighbor> neighbors, std::map<OctIndex, Decision> uncleDecisions);
  void interpolateAndSend(int);
  void interpolateAndSendToNephew(int, OctIndex);
  void refine();
  void interpolate(float*, std::vector<float>&, int, int, int, int, int, int);
  void refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, int zstart, int zend, float x_min, float y_min, float z_min);

  /*Load Balancing functions*/
  void startLdb();
  void ResumeFromSync();
  void UserSetLBLoad() { setObjTime((isRefined?0:1)); }

  bool isRoot();

  // AMR3D
  int amr3d_i;
  int ichild;
  void printData();
};

class MeshManager : public CBase_MeshManager {
  int workUnitCount;
  std::map<int, std::pair<float, float> > qdtimes;
  std::map<int, std::pair<float, float> > remeshtimes;
  std::map<int, int> workUnits;
  std::map<int, int> minLoad;
  std::map<int, int> maxLoad;
  std::map<int, float> avgLoad;

  double compute_time_sum;
  int compute_time_cnt;
  double decision_time_sum;
  int decision_time_cnt;

 public:
  float ****delu, ****delua;
  float delu2[numDims2], delu3[numDims2], delu4[numDims2];
#ifdef USE_GPU
  float* d_delu;
  float* d_delua;
#endif

  MeshManager_SDAG_CODE
  MeshManager();
  MeshManager(CkMigrateMessage *m);
  ~MeshManager();
  void pup(PUP::er &p){}
  void incrementWorkUnitCount(int);
  void recordQdTime(int iter, float a, float b){
    if (qdtimes.find(iter) == qdtimes.end())
      qdtimes[iter] = std::pair<float, float>(std::numeric_limits<float>::min(), std::numeric_limits<float>::max());
    qdtimes[iter].first = std::max(qdtimes[iter].first, a);
    qdtimes[iter].second = std::min(qdtimes[iter].second, b);
  }

  void recordRemeshTime(int iter, float a, float b){
    if (remeshtimes.find(iter) == remeshtimes.end())
      remeshtimes[iter] = std::pair<float, float>(0, std::numeric_limits<float>::max());
    remeshtimes[iter].first = std::max(remeshtimes[iter].first, a);
    remeshtimes[iter].second = std::min(remeshtimes[iter].second, b);
  }

  void processQdTimes(map<int, pair<float, float> > peQdtimes, map<int, pair<float, float> > peRemeshtimes, map<int, int> peWorkunits, map<int, int> peminLoad, map<int, int> pemaxLoad, map<int, float> peavgLoad);

  void printLogs();

  void reduceWorkUnits();
  void meshGenerationPhaseIsOver();

  void addComputeTime(double time) {
    compute_time_sum += time;
    compute_time_cnt++;
  }
  void addDecisionTime(double time) {
    decision_time_sum += time;
    decision_time_cnt++;
  }
};

extern CProxy_MeshManager mesh_manager;

#endif // ADVECTION_H_
