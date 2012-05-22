#ifdef LOGGER
#define VB(x) do { x } while(false)
#else
#define VB(x) do { } while(false)
#endif

#include "Constants.h"
#include "QuadIndex.h"
#include "pup_stl.h"
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
    else if(quad==3){
        if(dir==RIGHT)
            return RIGHT_DOWN;
        else if(dir==DOWN)
            return DOWN_RIGHT;
    }
}

inline void getChildren(QuadIndex myIndex, DIR dir, QuadIndex& q1, QuadIndex& q2){
  string str = myIndex.getIndexString();
  string q1str = str, q2str = str;

    if(dir==LEFT){
        q1str += "01";
        q2str += "10";
    }
    else if(dir==RIGHT){
        q1str += "00";
        q2str += "11";
    }
    else if(dir==UP){
        q1str += "01";
        q2str += "00";
    }
    else if(dir==DOWN){
        q1str += "10";
        q2str += "11";
    }

    q1 = QuadIndex(q1str.c_str());
    q2 = QuadIndex(q2str.c_str());
    return;
}


int inline wrap(int item, int max_size){
    return item%max_size;
}

inline const char* map_child(int child){
    if(child == LEFT_UP || child == UP_LEFT)
        return "01";
    else if(child == LEFT_DOWN || child==DOWN_LEFT)
        return "10";
    else if(child == DOWN_RIGHT || child==RIGHT_DOWN)
        return "11";
    else if(child == RIGHT_UP || child == UP_RIGHT)
        return "00";
    else return "-1";
}

class ChildDataMsg;

typedef TerminationChare<CProxy_Advection, QuadIndex> AdvTerm;

class Advection: public CBase_Advection,
  public AdvTerm {
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

        /*Phase1 DataStructures*/
	DECISION decision;
        DECISION nbr_decision[NUM_NEIGHBORS+2*NUM_NEIGHBORS];//Keeps the state of the neighbors
        DECISION child_decision[NUM_CHILDREN];

        bool hasReset;
        bool hasInitiatedPhase1;
        bool hasInitiatedPhase2;
	bool parentHasAlreadyMadeDecision;//to be used by a parent
	bool hasReceivedParentDecision;
        bool hasCommunicatedSTAY, hasCommunicatedREFINE;

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

        int iterations, cycle;
        bool isActive, shouldDestroy;
        
        double up;
        double un;
        double myt, mydt;
        double dx, dy, nx, ny;
        double xmin, xmax, ymin, ymax;

        void mem_allocate(double* &p, int size);
        void mem_allocate_all();
        QuadIndex getRefinedNeighbor(int NBR);
        int getSourceDirection(int NBR);
        double* getGhostBuffer(int dir);
        int getGhostCount(int dir);

        ~Advection();
        void free_memory(){/* Place Holder for calling Advection destructor - Advection::~Advection();*/}
        
        /*Constructors*/
        Advection(double, double, double, double);
        Advection(InitRefineMsg*);
  Advection() : AdvTerm(thisProxy, thisIndex, true) {advection(); 
          ckout << thisIndex.getIndexString().c_str() << " created 3" << endl;
  }
  Advection(CkMigrateMessage* m) : AdvTerm(thisProxy, thisIndex, true) {__sdag_init();
  ckout << thisIndex.getIndexString().c_str() << " created 4" << endl;
}
        
        void advection();// common function for initialization

        void printState();
        void done();
        void pup(PUP::er &p);
        
        /*Computation Methods*/
        /*void doStep(){
            begin_iteration();
        }*/
        /*void receiveGhosts(int iter, int dir, int width, double *u){
            imsg=0;
            process(iter, dir, width, u);
            sendReadyData();
        }*/
        void begin_iteration();
        void process(int, int, int, double*);
        void compute_and_iterate();
        void iterate();

        /*Phase1 Entry Methods*/
        DECISION getGranularityDecision();

        void doMeshRestructure();
        void resetMeshRestructureData();
        void communicatePhase1Msgs();
        void informParent(int, DECISION);
        void recvParentDecision();
        //void recvNeighborDecision(DIR);
        //void recvStatusUpdateFromParent(int);
        void exchangePhase1Msg(int, DECISION);
        void phase1DoneQD();

        /*Phase2 entry methods*/
        void setNbrStatus(int, ChildDataMsg*);
        void sendReadyData();
        void sendReadyData2RefiningNeighbors();
        // Returns whether a message was sent
        bool sendGhost(int,bool);
        QuadIndex lastSent;
        void doPhase2();

        void recvChildData(ChildDataMsg*);
        void interpolateAndSend(int);
        void refine();
        void interpolate(double*, vector<double>&, int, int, int, int);
        void refineChild(const char* sChild, int xstart, int xend, int ystart, int yend, double xmin, double ymin);
        template<class T>
        void print_Array(T*,int,int);

        /*LiveViz*/
        void requestNextFrame(liveVizRequestMsg*);

        void rootTerminated2();
};

class InitRefineMsg: public CMessage_InitRefineMsg{

    public:
        double dx, dy, myt, mydt, xmin, ymin, *refined_u;
        int iterations;
        bool *parent_nbr_exists;
        bool *parent_nbr_isRefined;
        DECISION *parent_nbr_decision;

        InitRefineMsg(){};
        InitRefineMsg(double dx, double dy, double myt, double mydt, double xmin, double ymin, int iterations, 
                        vector<double>& refined_u, bool *nbr_exists, bool *nbr_isRefined, DECISION *nbr_decision);
};

class ChildDataMsg: public CMessage_ChildDataMsg{
    public:
        int childNum;
        double myt, mydt; int iterations;
	double *child_u;
        bool *child_nbr_exists;
        bool *child_nbr_isRefined;
        DECISION *child_nbr_decision;
        
        ChildDataMsg(int cnum, double myt, double mydt, int iterations, double* u, bool* nbr_exists, bool* nbr_isRefined, DECISION* nbr_decision);
};

class PerProcessorChare : public CBase_PerProcessorChare{
    public:
        PerProcessorChare();
};
