#include "QuadIndex.h"
#include "Advection.decl.h"

class Advection: public CBase_Advection{
Advection_SDAG_CODE
    public:
        //tree information
        bool isdummy;
        bool hasRealChildren;
        bool hasDummyChildren;

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

        int iterations;
        
        double up;
        double un;

        void mem_allocate(double* &p, int size);
        void mem_allocate_all();
        
        Advection(bool, bool, bool);
        Advection(){advection();}
        Advection(CkMigrateMessage* m) {__sdag_init();}
        
        void advection();// common function for initialization
        void refine();
        void inform_nbr_of_refinement();

        void setDummy(){ isdummy = true;}
        void setReal(){ isdummy = false;}
        
        void pup(PUP::er &p);
        ~Advection();

        void begin_iteration();
        void process(int, int, double*);
        void compute_and_iterate();
        void iterate();
        void requestNextFrame(liveVizRequestMsg*);
};
