#include "Constants.h"
#include "QuadIndex.h"
#include "Advection.decl.h"
#include "pup_stl.h"

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
    else{
        if(quad==3)
            return LEFT_DOWN;
        else if(nbr == DOWN)
            return UP_RIGHT;
    }
    return -1;
}

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
    if(dir==LEFT){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString(), "01"));
        q2 = *new QuadIndex(strcat(myIndex.getIndexString(), "10"));
    }
    else if(dir==RIGHT){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString(), "00"));
        q2 = *new QuadIndex(strcat(myIndex.getIndexString(), "11"));
    }
    else if(dir==UP){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString(), "01"));
        q2 = *new QuadIndex(strcat(myIndex.getIndexString(), "00"));
    }
    else if(dir==DOWN){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString(), "10"));
        q2 = *new QuadIndex(strcat(myIndex.getIndexString(), "11"));
    }
    return;
}


int inline wrap(int item, int max_size){
    return item%max_size;
}

inline char* map_child(int child){
    if(child == LEFT_UP || child == UP_LEFT)
        return "01";
    else if(child == LEFT_DOWN || DOWN_LEFT)
        return "10";
    else if(child == DOWN_RIGHT || RIGHT_DOWN)
        return "11";
    else if(child == RIGHT_UP || child == UP_RIGHT)
        return "00";
    else return "-1";
}

class ChildDataMsg;

class Advection: public CBase_Advection{
Advection_SDAG_CODE
    public:
        //tree information
        bool exists;
        bool isRefined;
        bool isGrandParent;
        
        bool child_isRefined[NUM_CHILDREN];

        bool nbr_exists[NUM_NEIGHBORS];
        bool nbr_isRefined[NUM_NEIGHBORS];
        bool nbr_dataSent[NUM_NEIGHBORS];
        DECISION nbr_decision[NUM_NEIGHBORS+2*NUM_NEIGHBORS];//Keeps the state of the neighbors
        DECISION child_decision[NUM_CHILDREN];
        bool hasReceivedStatusFromParent[NUM_NEIGHBORS];

        set<int> hasReceived;
        bool hasInitiatedPhase1;
        bool hasInitiatedPhase2;
	bool parentHasAlreadyMadeDecision;
	bool hasReceivedParentDecision;

	bool hasAllocatedMemory;//for Use of Node who is going to derefine
	
        QuadIndex nbr[4], parent;
        int xc, yc;

        //data
        int imsg;

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
        
        double up;
        double un;
        double myt, mydt;
        double dx, dy;
        double xmin, xmax, ymin, ymax;

        void mem_allocate(double* &p, int size);
        void mem_allocate_all();

	DECISION decision;
        ~Advection();
        void free_memory(){/* Place Holder for calling Advection destructor - Advection::~Advection();*/}
        
        Advection(bool, bool, double, double, double, double);
        Advection(InitRefineMsg*);
        Advection(){advection();}
        Advection(CkMigrateMessage* m) {__sdag_init();}
        
        void advection();// common function for initialization

        void printState();
        void pup(PUP::er &p);

        void begin_iteration();
        void process(int, int, int, double*);
        void interpolateAndSend(int);
        void compute_and_iterate();
        void iterate();
        void refine();
        void interpolate(double*, double*, int, int, int, int);
        void requestNextFrame(liveVizRequestMsg*);

        void sendReadyData();
        void sendGhost(int);
        void doMeshRestructure();
        void communicateRefinement();
        void informParent(int, DECISION);
        void recvParentDecision();
        void recvNeighborDecision(DIR);
        void recvStatusUpdateFromParent(int);
        void exchangePhase1Msg(int);
        void recvChildData(ChildDataMsg*);
        DECISION getGranularityDecision();
        void doPhase2();

        void setNbrStatus(int, ChildDataMsg*);
};

class InitRefineMsg: public CMessage_InitRefineMsg{

    public:
        double dx, dy, myt, mydt, *refined_u;
        int iterations;
        bool parent_nbr_exists[NUM_NEIGHBORS];
        bool parent_nbr_isRefined[NUM_NEIGHBORS];
        bool parent_nbr_decision[3*NUM_NEIGHBORS];
        InitRefineMsg(double dx, double dy, double myt, double mydt, int iterations, double *refined_u, bool nbr_exists[NUM_NEIGHBORS], bool nbr_isRefined[NUM_NEIGHBORS], DECISION nbr_decision[3*NUM_NEIGHBORS]);
};

class ChildDataMsg: public CMessage_ChildDataMsg{
    public:
        int childNum;
	double *child_u;
        bool child_nbr_exists[NUM_NEIGHBORS];
        bool child_nbr_isRefined[NUM_NEIGHBORS];
        DECISION child_nbr_decision[NUM_NEIGHBORS];

        ChildDataMsg(){}
};
