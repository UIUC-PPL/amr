
class Main:public CBase_Main{

    public:
        int num_chares;
        int iterations;

        double io_outnum, io_tnext, io_tout; 

        Main(CkArgMsg* m);
        void printTreeInformation(CkVec<QuadIndex>);
        void terminate();
        void startMeshGeneration();
        void startRunning();
        void reportCascadeStats(int *cascade_lengths, int size);
        void qdlatency(double* elems, int size);
	void remeshlatency(double* elems, int size);
        void totalWorkUnits(int total);
};
