#include "Headers.h"
#include <assert.h>

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Advection qtree;

/* readonly */ int array_height, array_width, array_depth;
/* readonly */ int block_height, block_width, block_depth;
/* readonly */ int num_chare_rows, num_chare_cols, num_chare_Zs;
/* readonly */ int min_depth, max_depth;

/* readonly */ int max_iters, refine_freq, lb_freq;

/* readonly */ bool verbose;

/* readonly */ float x_min, x_max, y_min, y_max, z_min, z_max;
/* readonly */ float dx, dy, dz, vx, vy, vz;
/* readonly */ float apx, anx, apy, any, apz, anz;
/* readonly */ float x_ctr, y_ctr, z_ctr, radius;
/* readonly */ float tmax, t, dt, cfl;

/* readonly */ float start_time, end_time;

Main::Main(CkArgMsg* m) {
  mainProxy = thisProxy;

  // Set default parameters
  array_height = array_width = array_depth = 128;
  block_height = block_width = block_depth = 32;
  max_depth = 3;
  max_iters = 30;
  refine_freq = 3;
  lb_freq = 9;
  verbose = false;

  // Process command line arguments
  int c;
  while ((c = getopt(m->argc, m->argv, "a:b:d:i:r:l:vh")) != -1) {
    switch (c) {
      case 'a':
        array_height = array_width = array_depth = atoi(optarg);
        break;
      case 'b':
        block_height = block_width = block_depth = atoi(optarg);
        break;
      case 'd':
        max_depth = atoi(optarg);
        if (max_depth >= 11) {
          ckerr << "Depth is too large for bitvector index" << endl;
          CkExit();
        }
        break;
      case 'i':
        max_iters = atoi(optarg);
        break;
      case 'l':
        lb_freq = atoi(optarg);
        break;
      case 'r':
        refine_freq = atoi(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      case 'h':
        ckout << "\n[AMR Advection Mini-App Options]\n\n"
              << "\t-a\ttotal 3D array size\n"
              << "\t-b\tblock size per chare\n"
              << "\t-d\tmaximum depth of refinement\n"
              << "\t-i\tnumber of iterations\n"
              << "\t-r\trefinement frequency\n"
              << "\t-l\tload balancing frequency\n"
              << endl;
        CkExit();
      default:
        CkExit();
    }
  }

  for (int i = optind; i < m->argc; i++)
    printf("Non-option argument %s\n", m->argv[i]);

  // Check if array size is divisible by block size
  if (array_width < block_width || array_width % block_width != 0) {
    ckout << "Array size (" << array_width << ") should be divisible by block size ("
          << block_width << ")" << endl;
    CkExit();
  }

  // Check if load balancing frequency is a multiple of refine frequency
  if (lb_freq < refine_freq || lb_freq % refine_freq != 0) {
    ckerr << "Load balancing frequency (" << lb_freq << ") should be a multiple of"
          << "refine frequency (" << refine_freq << ")" << endl;
    CkExit();
  }

  // Set number of chares per dimension
  num_chare_rows = num_chare_cols = num_chare_Zs = array_width / block_width;
  int num_chares = num_chare_rows * num_chare_cols * num_chare_Zs;

  // Set minimum depth
  float fdepth = log(num_chares) / log(NUM_CHILDREN);
  min_depth = (fabs(fdepth - ceil(fdepth)) < 0.000001) ? ceil(fdepth) : floor(fdepth);
  if (min_depth > max_depth) {
    ckerr << "Minimum depth > maximum depth: try increasing maximum depth" << endl;
    CkExit();
  }
  if (min_depth == 0) {
    ckerr << "Minimum depth is 0, 1 or more chares are required to run the simulation" << endl;
    CkExit();
  }

  // Initialize constants
  x_min = 0; x_max = 1;
  y_min = 0; y_max = 1;
  z_min = 0; z_max = 1;

  dx = (x_max - x_min) / float(array_width);
  dy = (y_max - y_min) / float(array_height);
  dz = (z_max - z_min) / float(array_depth);

  vx = 0.0; vy = 0.0; vz = 0.1;

  apx = std::max(vx, (float)0.0);
  anx = std::min(vx, (float)0.0);
  apy = std::max(vy, (float)0.0);
  any = std::min(vy, (float)0.0);
  apz = std::max(vz, (float)0.0);
  anz = std::min(vz, (float)0.0);

  x_ctr = 0.5; y_ctr = 0.5; z_ctr = 0.5;
  radius = 0.2;

  t = 0; tmax = 10000;
  cfl = 0.4;
  dt = min(min(dx, dy), dz) / sqrt(vx * vx + vy * vy + vz * vz) * cfl;
  dt /= pow(2., max_depth - min_depth);
  if ((t + dt) >= tmax)
    dt = tmax - t;
  t = t + dt;

  CkPrintf("\n===== Welcome to the Charm++ AMR Mini-App =====\n"
           "* Array dimension: %d x %d x %d\n"
           "* Block dimension: %d x %d x %d\n"
           "* Initial chares: %d x %d x %d\n"
           "* Minimum depth: %d\n"
           "* Maximum depth: %d\n"
           "* Number of iterations: %d\n"
           "* Refinement frequency: %d\n"
           "* Load balancing frequency: %d\n\n",
           array_width, array_height, array_depth,
           block_width, block_height, block_depth,
           num_chare_rows, num_chare_cols, num_chare_Zs,
           min_depth, max_depth, max_iters, refine_freq, lb_freq);

  // Create tree of chares
  CProxy_AdvMap map = CProxy_AdvMap::ckNew();
  CkArrayOptions opts;
  opts.setMap(map);
  qtree = CProxy_Advection::ckNew(opts);
  for (int i = 0; i < num_chare_rows; ++i) {
    for (int j = 0; j < num_chare_cols; ++j) {
      for (int k = 0; k < num_chare_Zs; ++k) {
        qtree[OctIndex(i, j, k, min_depth)].insert(x_min, x_max, y_min, y_max, z_min, z_max);
      }
    }
  }
  qtree.doneInserting();

  // Begin simulation
  CkStartQD(CkCallback(CkIndex_Main::startMeshGeneration(), thisProxy));
  mesh_manager = CProxy_MeshManager::ckNew();
}

void Main::startMeshGeneration() {
  start_time = CkWallTimer();
  qtree.iterate();
}

void Main::terminate() {
  ckout << "\nSimulation time: " << CkWallTimer() - start_time << " s" << endl;
  mesh_manager.reduceWorkUnits();
}

void Main::totalWorkUnits(int total) {
  CkPrintf("Total work units: %d\n", total);
  mesh_manager.reduceQdTimes();
}

// Map chares to PEs
struct AdvMap : public CBase_AdvMap {
  int bits;
  AdvMap() : bits(log2(CkNumPes())) { }

  void pup(PUP::er &p) { bits = log2(CkNumPes()); }
  AdvMap(CkMigrateMessage* m): CBase_AdvMap(m), bits(log2(CkNumPes())) {}

  int procNum(int arrayHdl, const CkArrayIndex& i) {
    const unsigned long long golden_ratio_prime_64 = 0x9e37fffffffc0001ULL;
    int numPes = CkNumPes();
    const OctIndex& idx = *reinterpret_cast<const OctIndex*>(i.data());
    int baseBits = 8;

    unsigned long long val = idx.bitVector >> (sizeof(unsigned int)*8 - baseBits);
    unsigned long long hash = golden_ratio_prime_64 * val;

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
