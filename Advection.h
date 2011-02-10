#include "Constants.h"
#include "QuadIndex.h"
#include "Advection.decl.h"

class Advection: public CBase_Advection{
Advection_SDAG_CODE
    public:
        //tree information
        bool exits;
        bool isRefined;
        
        bool nbr_exists[NUM_NEIGHBORS];
        bool nbr_isRefined[NUM_NEIGHBORS];

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
        void mem_allocate(double* &p, int size);
        void mem_allocate_all();
        ~Advection();
        void free_memory(){/* Place Holder for calling Advection destructor - Advection::~Advection();*/}
        
        Advection(bool, bool, bool);
        Advection(){advection();}
        Advection(CkMigrateMessage* m) {__sdag_init();}
        
        void advection();// common function for initialization
        void refine();
        void derefine();
        void inform_nbr_of_refinement(int);
        void inform_nbr_hasDummyChildren(int inbr);
        void destroyChildren();
        void inform_nbr_hasNoChildren(int inbr);
        void inform_child_hasDummyChildren(int cindx){
            child_hasDummyChildren[cindx] = true;
        }

        void setDummy();
        void setReal();
        void manage_memory_RealToDummy();
        void manage_memory_DummyToReal();
        
        void printState();
        void pup(PUP::er &p);

        void begin_iteration();
        void process(int, int, int, double*);
        void compute_and_iterate();
        void iterate();
        void requestNextFrame(liveVizRequestMsg*);
};
