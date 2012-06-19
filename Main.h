
class Main:public CBase_Main{

    public:
        int num_chares;
        int iterations;

        double io_outnum, io_tnext, io_tout; 

        Main(CkArgMsg* m);
        void printTreeInformation(CkVec<QuadIndex>);
        void terminate();
        void startRunning();
};
