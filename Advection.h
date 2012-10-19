#if !defined(ADVECTION_H)
#define ADVECTION_H

#ifdef LOGGER
#define VB(x) do { x } while(false)
#else
#define VB(x) do { } while(false)
#endif

#include "Advection.decl.h"

int inline myDirectionWrtUncle(int quad, int nbr){
    switch(quad){
        case 0: return (nbr==RIGHT)?LEFT_UP:DOWN_RIGHT; break;
        case 1: return (nbr==LEFT)?RIGHT_UP:DOWN_LEFT;  break;
        case 2: return (nbr==LEFT)?RIGHT_DOWN:UP_LEFT; break;
        case 3: return (nbr==RIGHT)?LEFT_DOWN:UP_RIGHT; break;
        default: CkAbort("invalid quad#");
    };
}

//returns the direction of the neighbor in direction 'dir' with repsect to the parent
int inline nbrDirectionWrtParent(int quad, int dir){
  switch(quad){
      case 0: return (dir==RIGHT)?RIGHT_UP:UP_RIGHT;    break;
      case 1: return (dir==LEFT)?LEFT_UP:UP_LEFT;       break;
      case 2: return (dir==LEFT)?LEFT_DOWN:DOWN_LEFT;   break;
      case 3: return (dir==RIGHT)?RIGHT_DOWN:DOWN_RIGHT;break;
      default: CkAbort("invalid quad num");
  };
}

inline void getChildrenInDir(QuadIndex myIndex, int dir, QuadIndex& q1, QuadIndex& q2){
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

inline int childDir2Quadrant(int childDir){
  switch(childDir){
    case RIGHT_UP:      case UP_RIGHT:      return 0;
    case LEFT_UP:       case UP_LEFT:       return 1;
    case LEFT_DOWN:     case DOWN_LEFT:     return 2;
    case DOWN_RIGHT:    case RIGHT_DOWN:    return 3;
  }
}

//class ChildDataMsg;

/*typedef TerminationChare<CProxy_Advection, QuadIndex> AdvTerm;*/

class Advection: public CBase_Advection/*, public AdvTerm */{
  Advection_SDAG_CODE
    public:
        
  std::ofstream logFile;
  std::ofstream outFile;
  //tree information
  bool isRefined;
  int depth;

  bool child_isRefined[NUM_CHILDREN];
  bool isGrandParent();

  bool nbr_exists[NUM_NEIGHBORS];
  bool nbr_isRefined[NUM_NEIGHBORS];
  bool nbr_dataSent[3*NUM_NEIGHBORS];
  std::set<int> hasReceived;
  bool hasReceivedFromDir(int dir);
  bool hasReceivedFromAroundCorner(int aroundCorner);

  /*Phase1 DataStructures*/
  DECISION decision;
  DECISION nbr_decision[NUM_NEIGHBORS+2*NUM_NEIGHBORS];//Keeps the state of the neighbors
  DECISION child_decision[NUM_CHILDREN];

  bool parentHasAlreadyMadeDecision;//to be used by a parent
  bool hasReceivedParentDecision;

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
  Advection(double dx, double dy,
            double myt, double mydt, double xmin, double ymin,
            int meshGenIterations_, int iterations_, 
            std::vector<double> refined_u, bool *parent_nbr_exists, 
            bool *parent_nbr_isRefined, DECISION *parent_nbr_decision);

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
  void pup(PUP::er &p);

  /* initial mesh generation*/       
  //void updateMesh();
  void applyInitialCondition();
        
  void process(int, int, int, double*);
  void compute();
  void iterate();

  /*Phase1 Entry Methods*/
  void makeGranularityDecisionAndCommunicate();
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
  void setNbrStateUponCoarsening(int, bool*, bool*, DECISION*);
  void sendReadyData();
  // Returns whether a message was sent
  void sendGhost(int);
  void doPhase2();
  void updateMeshState();

  void recvChildData(int, double, double, int, int, std::vector<double>, bool*, bool*, DECISION*);
  void interpolateAndSend(int);
  void refine();
  void interpolate(double*, std::vector<double>&, int, int, int, int);
  void refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, double xmin, double ymin);

  /*Load Balancing functions*/
  void startLdb();
  void ResumeFromSync();
  void UserSetLBLoad();

  void rootTerminated2();

  bool isRoot();
};

/*class InitRefineMsg: public CMessage_InitRefineMsg{

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
};*/

/*class ChildDataMsg: public CMessage_ChildDataMsg{
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
};*/

class AdvectionGroup : public CBase_AdvectionGroup {
  std::vector<int> cascades;
  std::vector<double> qdlatencies, remeshlatencies;
  int workUnitCount;
  int numNotificationsExpected, numNotificationsRecvd;
  int nChanges;
  public:
  bool meshUpdated;

 public:
  AdvectionGroup();

  void incrementWorkUnitCount();
  void recordCascade(int iteration, int length);
  void collectCascades(CkCallback cb);
  void recordQDLatency(int iteration, double latency);
  void recordRemeshLatency(int iteration, double latency);
  void reduceWorkUnits();
  void reduceLatencies();
  void meshGenerationPhaseIsOver();
  void notifyMeshUpdate(DECISION);
  void meshUpdateReductionClient(int);
  void resetMeshUpdateCounters();
};

extern CProxy_AdvectionGroup ppc;

#endif // ADVECTION_H
