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
        
        bool nbr_exists[NUM_NEIGHBORS];
        bool nbr_isRefined[NUM_NEIGHBORS];
        bool nbr_dataSent[NUM_NEIGHBORS];
        set<int> hasReceived;

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
        ~Advection();
        void free_memory(){/* Place Holder for calling Advection destructor - Advection::~Advection();*/}
        
        Advection(bool, bool, double, double, double, double);
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
        void requestNextFrame(liveVizRequestMsg*);
};
