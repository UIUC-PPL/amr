#include "Headers.h"

using std::max;
using std::min;

CProxy_Advection qtree;
CProxy_Main mainProxy;

int array_height;
int array_width;

int num_chare_rows;
int num_chare_cols;

int block_height;
int block_width;

int min_depth, max_depth;
int max_iterations, refine_frequency;

const char amrRevision[] = INQUOTES(AMR_REVISION);

double xmin, xmax, ymin, ymax;
double xctr, yctr, radius;
double dx, dy, v;
double ap, an;
double tmax, t, dt, cfl;
bool inInitialMeshGenerationPhase;

double start_time, end_time;

Main::Main(CkArgMsg* m){
  ckout<<"Running amr code revision: "<<amrRevision<<endl;

  mainProxy = thisProxy;
  iterations = 0;

  if(m->argc < 4){
    ckout << "Usage: " << m->argv[0] << "[max_depth] [block_size] [iterations] [array_dim]?" << endl; 
    CkExit();
  }

  if (m->argc >= 5) {
    array_height = array_width = atoi(m->argv[4]);
  } else {
    array_height = array_width = 256;
  }
  
  block_height = block_width = atoi(m->argv[2]);
  max_iterations = atoi(m->argv[3]);

  if(array_width%block_width < 0 || array_width < block_width){
    ckout << "Incompatible arguments" << endl;
    CkExit();
  }

  num_chare_rows = num_chare_cols = array_height/block_height;
    
  /******** Do Some Initialization *********/
  xmin = 0;
  xmax = 1;
  ymin = 0;
  ymax = 1;
  t= 0;
  tmax = 10000;
  cfl = 0.4;
  v = 0.1;
  refine_frequency = 3;
  inInitialMeshGenerationPhase = true;

  dx = (xmax - xmin)/double(array_width);
  dy = (ymax - ymin)/double(array_height);
  xctr = 0.3;
  yctr = 0.5;
  radius = 0.2;

  ap = max(v, 0.0);
  an = min(v, 0.0);

  dt = min(dx,dy)/v * cfl;
  dt /= pow(2., max_depth - min_depth);
  if ((t + dt) >= tmax )
    dt = tmax - t;
  t = t+dt;

  /*****End Initialization **********/
  CProxy_AdvMap map = CProxy_AdvMap::ckNew();
  CkArrayOptions opts;
  opts.setMap(map);
  qtree = CProxy_Advection::ckNew(opts);

  //save the total number of worker chares we have in this simulation
  num_chares = num_chare_rows*num_chare_cols;
  double fdepth = (log(num_chares)/log(4));
  int depth = (fabs(fdepth - ceil(fdepth)) < 0.000001)?ceil(fdepth):floor(fdepth);
  min_depth = depth;
  CkAssert(min_depth >= 4);

  // To maintain the semantics of "max_depth" that set it relative to
  // a grid fo 256, offset by 4
  max_depth = atoi(m->argv[1]) + min_depth - 4;

  CkPrintf("Running Advection on %d processors with (%d,%d) elements, maxDepth = %d, blockSize = %d, maxIter = %d\n",
           CkNumPes(), array_width, array_height, max_depth, block_height, max_iterations);

  for (int i = 0; i < num_chare_rows; ++i)
    for (int j = 0; j < num_chare_cols; ++j)
      qtree[QuadIndex(i, j, min_depth)].insert(xmin, xmax, ymin, ymax);

  qtree.doneInserting();

  CkStartQD(CkCallback(CkIndex_Main::startMeshGeneration(), thisProxy));
  ppc = CProxy_AdvectionGroup::ckNew();

}

void Main::startMeshGeneration() {
  start_time = CkWallTimer();
  qtree.iterate();
}

void Main::terminate(){
  ckout << "simulation time: " << CkWallTimer() - start_time << " s" << endl;
  ppc.reduceWorkUnits();
}

void Main::totalWorkUnits(int total) {
  CkPrintf("total work units = %d\n", total);
  CkExit();
}

#define GOLDEN_RATIO_PRIME_64 0x9e37fffffffc0001ULL

struct AdvMap : public CBase_AdvMap {
  int bits;
  AdvMap() : bits(log2(CkNumPes())) { }

  int procNum(int arrayHdl, const CkArrayIndex& i) {
    int numPes = CkNumPes();
    const QuadIndex& idx = *reinterpret_cast<const QuadIndex*>(i.data());
    int baseBits = 8;

    unsigned long long val = idx.bitVector >> (sizeof(unsigned int)*8 - baseBits);
    unsigned long long hash = GOLDEN_RATIO_PRIME_64 * val;

    int basePE = hash >> (64 - bits);

    unsigned long validBits = idx.bitVector & ((1L << 24) - 1);
    validBits += (1L << 22);
    unsigned long offset = validBits >> (sizeof(unsigned int)*8 - idx.nbits);
    offset += (idx.nbits == 8);

    int pe = (basePE + offset - 1) % numPes;

    return pe;
  }
};

#include "Main.def.h"
