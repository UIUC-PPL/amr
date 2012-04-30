#include "boost/assign.hpp"
#include "Advection.h"
#include "Main.decl.h"
#include <vector>

class Main:public CBase_Main{

    public:
        CProxy_Advection qtree;
        int num_chares;
        int iterations;

        double io_outnum, io_tnext, io_tout; 

        Main(CkArgMsg* m);
        void printTreeInformation(CkVec<QuadIndex>);
        void initUtils();
        void terminate();
        void startRunning();
        void kickPhase2();
};
