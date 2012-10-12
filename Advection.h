
#if !defined(ADVECTION_H)
#define ADVECTION_H

#ifdef LOGGER
#define VB(x) do { x } while(false)
#else
#define VB(x) do { } while(false)
#endif

#include "Advection.decl.h"

int inline map_nbr(int quad, int nbr){
  if(quad==0){
    if(nbr == RIGHT)
      return LEFT_UP;
    else if(nbr == UP)
      return DOWN_RIGHT;
  }
  else if(quad==1){
    if(nbr == LEFT)
      return RIGHT_UP;
    else if(nbr == UP)
      return DOWN_LEFT;
  }
  else if(quad==2){
    if(nbr == LEFT)
      return RIGHT_DOWN;
    else if(nbr == DOWN)
      return UP_LEFT;
  }
  else if(quad==3){
    if(nbr==RIGHT)
      return LEFT_DOWN;
    else if(nbr == DOWN)
      return UP_RIGHT;
  }
  return -1;
}

//returns the direction of the neighbor in direction 'dir' with repsect to the parent
int inline getNbrDir(int quad, int dir){
  if(quad==0){
    if(dir==RIGHT)
      return RIGHT_UP;
    else if(dir==UP)
      return UP_RIGHT;
  }
  else if(quad==1){
    if(dir==LEFT)
      return LEFT_UP;
    else if(dir==UP)
      return UP_LEFT;
  }
  else if(quad==2){
    if(dir==LEFT)
      return LEFT_DOWN;
    else if(dir==DOWN)
      return DOWN_LEFT;
  }
  else{// if(quad==3){
    if(dir==RIGHT)
      return RIGHT_DOWN;
    else if(dir==DOWN)
      return DOWN_RIGHT;
  }
}

inline void getChildren(QuadIndex myIndex, DIR dir, QuadIndex& q1, QuadIndex& q2){
  int q1c, q2c;
  switch (dir) {
  case LEFT:  q1c = 1; q2c = 2; break;
  case RIGHT: q1c = 0; q2c = 3; break;
  case UP:    q1c = 1; q2c = 0; break;
  case DOWN:  q1c = 2; q2c = 3; break;
  }

  q1 = myIndex.getChild(q1c);
  q2 = myIndex.getChild(q2c);
  return;
}


int inline wrap(int item, int max_size){
  return item%max_size;
}

inline int map_child(int child){
  if(child == LEFT_UP || child == UP_LEFT)
    return 1;
  else if(child == LEFT_DOWN || child==DOWN_LEFT)
    return 2;
  else if(child == DOWN_RIGHT || child==RIGHT_DOWN)
    return 3;
  else if(child == RIGHT_UP || child == UP_RIGHT)
    return 0;
  else return -1;
}

class ChildDataMsg;

/*typedef TerminationChare<CProxy_Advection, QuadIndex> AdvTerm;*/

class Advection: public CBase_Advection/*, public AdvTerm */{
  Advection_SDAG_CODE
    public:
        
  ofstream logFile;
  ofstream outFile;
  //tree information
  bool isRefined;
  bool isGrandParent;
  int depth;

  bool child_isRefined[NUM_CHILDREN];

  bool nbr_exists[NUM_NEIGHBORS];
  bool nbr_isRefined[NUM_NEIGHBORS];
  bool nbr_dataSent[3*NUM_NEIGHBORS];
  set<int> hasReceived;
  bool hasReceivedFromDir(int dir);
  bool hasReceivedFromAroundCorner(int aroundCorner);

  /*Phase1 DataStructures*/
  DECISION decision;
  DECISION nbr_decision[NUM_NEIGHBORS+2*NUM_NEIGHBORS];//Keeps the state of the neighbors
  DECISION child_decision[NUM_CHILDREN];

  bool hasReset;
  bool parentHasAlreadyMadeDecision;//to be used by a parent
  bool hasReceivedParentDecision;

  bool hasAllocatedMemory;//for Use of Node who is going to derefine
  bool has_terminated; //used to inform parent about termination so that it can contribute to the final reduction

  QuadIndex nbr[4], parent;
  int xc, yc;

  //data
  double imsg;

  double* u;
  double* u2;
  double* u3;
  double *x;
  double *y;

  double *left_edge;
  double *right_edge;
        
  /* Required In Case of Dummy Nodes */
  double *top_edge;
  double *bottom_edge;

  int iterations;
  int meshGenIterations;
  bool shouldDestroy;
        
  double up;
  double un;
  double myt, mydt;
  double dx, dy, nx, ny;
  double xmin, xmax, ymin, ymax;
        
  double itBeginTime, remeshStartTime;
  void mem_allocate(double* &p, int size);
  void mem_allocate_all();
  QuadIndex getRefinedNeighbor(int NBR);
  int getSourceDirection(int NBR);
  double* getGhostBuffer(int dir);
  int getGhostCount(int dir);
  double lastIdleTimeQD;

  ~Advection();
  void free_memory(){/* Place Holder for calling Advection destructor - Advection::~Advection();*/}
        
  /*Constructors*/
  Advection(double, double, double, double);
  Advection(InitRefineMsg*);
  Advection() /*: AdvTerm(thisProxy, thisIndex, true) */{
    usesAutoMeasure = CmiFalse;
    //advection(); 
    //ckout << thisIndex.getIndexString().c_str() << " created 3" << endl;
  }
  Advection(CkMigrateMessage* m) /*: AdvTerm(thisProxy, thisIndex, true) */{
    usesAutoMeasure = CmiFalse;
    __sdag_init();
    //ckout << thisIndex.getIndexString().c_str() << " created 4" << endl;
  }
        
  void advection();// common function for initialization

  void printState();
  void done();
  void pup(PUP::er &p);

  /* initial mesh generation*/       
  void updateMesh();
  void applyInitialCondition();
        
  void process(int, int, int, double*);
  void compute_and_iterate();
  void iterate();

  /*Phase1 Entry Methods*/
  void doRemeshing();
  DECISION getGranularityDecision();

  //void doMeshRestructure();
  void resetMeshRestructureData();
  void updateDecisionState(int cascade_length, DECISION newDecision);
  void informParent(int, DECISION, int cascade_length);
  void recvParentDecision(int cascade_length);
  //void recvNeighborDecision(DIR);
  //void recvStatusUpdateFromParent(int);
  void exchangePhase1Msg(int, DECISION, int cascade_length);

  /*Phase2 entry methods*/
  void setNbrStatus(int, ChildDataMsg*);
  void sendReadyData();
  void sendReadyData2RefiningNeighbors();
  // Returns whether a message was sent
  bool sendGhost(int,bool);
  QuadIndex lastSent;
  void doPhase2();
  void updateNbrStatus();

  void recvChildData(ChildDataMsg*);
  void interpolateAndSend(int);
  void refine();
  void interpolate(double*, vector<double>&, int, int, int, int);
  void refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, double xmin, double ymin);
  template<class T>
    void print_Array(T*,int,int);

  /*Load Balancing functions*/
  void startLdb();
  void ResumeFromSync();
  void UserSetLBLoad();

  /*LiveViz*/
  void requestNextFrame(liveVizRequestMsg*);

  void rootTerminated2();

  bool isRoot();
};

class InitRefineMsg: public CMessage_InitRefineMsg{

 public:
  bool isInMeshGenerationPhase;
  double dx, dy, myt, mydt, xmin, ymin, *refined_u;
  int meshGenIterations, iterations;
  bool *parent_nbr_exists;
  bool *parent_nbr_isRefined;
  DECISION *parent_nbr_decision;

  InitRefineMsg(){};
  InitRefineMsg(bool isInMeshGenerationPhase, double dx, double dy, double myt, double mydt, double xmin, double ymin, 
                int meshGenIterations, int iterations, vector<double>& refined_u, bool *nbr_exists, bool *nbr_isRefined, DECISION *nbr_decision);
};

class ChildDataMsg: public CMessage_ChildDataMsg{
 public:
  bool isInMeshGenerationPhase;
  int childNum;
  double myt, mydt; 
  int meshGenIterations, iterations;
  double *child_u;
  bool *child_nbr_exists;
  bool *child_nbr_isRefined;
  DECISION *child_nbr_decision;
        
  ChildDataMsg(bool isInMeshGenerationPhase, int cnum, double myt, double mydt, int meshGenIterations, int iterations, double* u, bool* nbr_exists, bool* nbr_isRefined, DECISION* nbr_decision);
};

class PerProcessorChare : public CBase_PerProcessorChare {
  vector<int> cascades;
  std::vector<double> qdlatencies, remeshlatencies;
  int workUnitCount;

 public:
  PerProcessorChare();

  void incrementWorkUnitCount();
  void recordCascade(int iteration, int length);
  void collectCascades(CkCallback cb);
  void recordQDLatency(int iteration, double latency);
  void recordRemeshLatency(int iteration, double latency);
  void reduceWorkUnits();
  void reduceLatencies();
};

extern CProxy_PerProcessorChare ppc;

#endif // ADVECTION_H
