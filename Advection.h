#include "Constants.h"
#include "QuadIndex.h"
#include "Advection.decl.h"
#include "pup_stl.h"
#include "Messages.h"

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

inline void getChildren(QuadIndex myIndex, DIR dir, QuadIndex& q1, QuadIndex& q2){
    if(dir==LEFT){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString()), "01");
        q2 = *new QuadIndex(strcat(myIndex.getIndexString()), "10");
    }
    else if(dir==RIGHT){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString()), "00");
        q2 = *new QuadIndex(strcat(myIndex.getIndexString()), "11");
    }
    else if(dir==UP){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString()), "01");
        q2 = *new QuadIndex(strcat(myIndex.getIndexString()), "00");
    }
    else if(dir==DOWN){
        q1 = *new QuadIndex(strcat(myIndex.getIndexString()), "10");
        q2 = *new QuadIndex(strcat(myIndex.getIndexString()), "11");
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
        DIR nbr_decision[NUM_NEIGHBORS+2*NUM_NEIGHBORS];//Keeps the state of the neighbors
        DIR child_decision[NUM_CHILDREN];

        set<int> hasReceived;
        bool hasInitiatedPhase2;
	bool parentHasAlreadyMadeDecision;
	bool hasReceivedParentDecision;

	bool hasAllocatedMemory;

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
        void doPhase2();
        void requestNextFrame(liveVizRequestMsg*);
};

class InitRefineMsg: public CBase_InitRefineMsg{
    public:
        double dx, dy, myt, mydt, *refined_u;

};

class ChildDataMsg: public CBase_ChildDataMsg{
    public:
        int childNum;
	double *child_u;
	
	ChildDataMsg(int num, double *u){
	    childNum = num;
	    memcpy();
	}
};
