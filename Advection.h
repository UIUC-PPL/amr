#if !defined(ADVECTION_H)
#define ADVECTION_H

#include "Constants.h"

#include "Advection.decl.h"

class Neighbor {

  int upper_bound, lower_bound;
  bool refined;
  bool dataSent;
  Decision decision;
  Decision childDecisions[NUM_CHILDREN];
  int childBounds[2*NUM_CHILDREN];
  int dir;

public:

  Neighbor() : refined(false), dataSent(false), decision(INV),
               upper_bound(100), lower_bound(-1), dir(-1) {}

  Neighbor(int dir) : refined(false), dataSent(false), decision(INV),
                      upper_bound(100), lower_bound(-1), dir(dir) {}

  int getDir() { return dir; }

  bool isRefined() {
    return refined;
  }

  bool isDataSent() {
    return dataSent;
  }

  bool isConverged(int child = -1) {
    if (child == -1) {
      return lower_bound == upper_bound;
    } else {
      return childBounds[2*child+0] == childBounds[2*child+1];
    }
  }

  void setDataSent(bool value) {
    dataSent = value;
  }

  void setRefined(bool value) {
    refined = value;
  }

  void setBounds(int lower, int upper, int child = -1) {
    if (child == -1) {
      lower_bound = max(lower_bound, lower);
      upper_bound = min(upper_bound, upper);
    } else {
      childBounds[2*child+0] = max(childBounds[2*child+0], lower);
      childBounds[2*child+1] = min(childBounds[2*child+1], upper);
    }
  }

  int getUpperBound(int child = -1) {
    if (child == -1) {
      return upper_bound;
    } else {
      return childBounds[2*child+1];
    }
  }

  int getLowerBound(int child = -1) {
    if (child == -1) {
      return lower_bound;
    } else {
      return childBounds[2*child];
    }
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
    p|upper_bound;
    p|lower_bound;
    p|refined;
    p|dataSent;
    p|decision;
    p|dir;
    PUParray(p, childDecisions, NUM_CHILDREN);
    PUParray(p, childBounds,  2*NUM_CHILDREN);
  }

  void resetDecision() {
    decision = COARSEN;
    lower_bound =  -1;
    upper_bound = 100;
    for (int i = 0; i < NUM_CHILDREN; ++i) {
      childDecisions[i] = COARSEN;
      childBounds[2*i+0] =  -1;
      childBounds[2*i+1] = 100;
    }
  }
};

//class ChildDataMsg;

/*typedef TerminationChare<CProxy_Advection, OctIndex> AdvTerm;*/

class Advection: public CBase_Advection/*, public AdvTerm */{
  Advection_SDAG_CODE
    public:
        
  //tree information
  bool isRefined;
  bool isMaxRefined;
  int depth;

  bool child_isRefined[NUM_CHILDREN];
  bool isGrandParent();

  std::set<OctIndex> ghostReceived;

  /*Phase1 DataStructures*/
  Decision decision;
  Decision child_decision[NUM_CHILDREN];

  bool parentHasAlreadyMadeDecision;//to be used by a parent
  bool hasReceivedParentDecision;
  char fname[100];
  OctIndex  parent;
  std::map<OctIndex, Neighbor> neighbors;
  std::map<OctIndex, Decision> uncleDecisions;
  int xc, yc, zc;

  int lower_bound, upper_bound;

  //data
  float imsg;

  float* u;
  float* u2;
  float* u3;
  float *x;
  float *y;
  float *z;

  float *left_surface;
  float *right_surface;
  float *top_surface;
  float *bottom_surface;
  float *forward_surface;
  float *backward_surface;

  int iterations;
  int meshGenIterations;

  float up;
  float un;
  float myt, mydt;
  float dx, dy, dz, nx, ny, nz;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  
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

  ~Advection();

  /* Constructors */
  Advection() { usesAtSync = true; usesAutoMeasure = false; }
  Advection(CkMigrateMessage* m) : CBase_Advection(m){ 
    usesAutoMeasure = false;
    usesAtSync = true;
    VB(logfile.open(string("log/"+thisIndex.getIndexString()+"log").c_str(), std::ofstream::app););
    VB(logfile << "migrated to a new processor" << std::endl;)
  }
  Advection(float, float, float, float, float, float);
  Advection(float dx, float dy, float dz,
            float myt, float mydt,
            float xmin, float ymin, float zmin,
            int meshGenIterations_, int iterations_,
            std::vector<float> refined_u,
            std::map<OctIndex, Neighbor> neighbors);

  void initializeRestofTheData(); // common function initialization for
  
  void pup(PUP::er &p);

  void updateBounds(int new_lower_bound, int new_upper_bound, bool notify_neighbors = true);
  void notifyAllNeighbors(int cascade_length);

  /* initial mesh generation*/
  void applyInitialCondition();
  void process(int, int, int, int, float*);
  void compute();

  /*Phase1 Entry Methods*/
  void makeGranularityDecisionAndCommunicate();
  Decision getGranularityDecision();

  void resetMeshRestructureData();
  void prepareData4Exchange();
  void processPhase1Msg(int, int, Decision, int, int, int, OctIndex);

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
  void refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, int zstart, int zend, float xmin, float ymin, float zmin);

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

class AdvectionGroup : public CBase_AdvectionGroup {
  int workUnitCount;
  std::map<int, std::pair<float, float> > qdtimes;
  std::map<int, std::pair<float, float> > remeshtimes;
  std::map<int, int> workUnits;
  std::map<int, int> minLoad;
  std::map<int, int> maxLoad;
  std::map<int, float> avgLoad;
 public:
  AdvectionGroup_SDAG_CODE
  AdvectionGroup();
  AdvectionGroup(CkMigrateMessage *m);
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
};

extern CProxy_AdvectionGroup ppc;

#endif // ADVECTION_H
