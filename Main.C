#include "Headers.h"
#include <assert.h>

/* readonly */ CProxy_Main main_proxy;
/* readonly */ CProxy_MeshBlock mesh;
/* readonly */ int grid_y, grid_x, grid_z;
/* readonly */ int block_y, block_x, block_z;
/* readonly */ int n_chares_x, n_chares_y, n_chares_z;
/* readonly */ int min_depth, max_depth;
/* readonly */ int max_iters, refine_freq, lb_freq;
/* readonly */ bool verbose;
/* readonly */ float x_min, x_max, y_min, y_max, z_min, z_max;
/* readonly */ float dx, dy, dz, vx, vy, vz;
/* readonly */ float apx, anx, apy, any, apz, anz;
/* readonly */ float x_ctr, y_ctr, z_ctr, radius;
/* readonly */ float tmax, t, dt, cfl;

Main::Main(CkArgMsg* m) {
  main_proxy = thisProxy;

  // Set default parameters
  grid_x = grid_y = grid_z = 128;
  block_x = block_y = block_z = 32;
  max_depth = 3;
  max_iters = 30;
  refine_freq = 3;
  lb_freq = 9;
  verbose = false;

  // Process command line arguments
  int c;
  while ((c = getopt(m->argc, m->argv, "a:b:d:i:r:l:v")) != -1) {
    switch (c) {
      case 'a':
        grid_x = grid_y = grid_z = atoi(optarg);
        break;
      case 'b':
        block_x = block_y = block_z = atoi(optarg);
        break;
      case 'd':
        max_depth = atoi(optarg);
        if (max_depth >= 11) {
          CkAbort("Depth is too large for bitvector index");
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
      default:
        CkExit();
    }
  }

  for (int i = optind; i < m->argc; i++) {
    CkPrintf("Non-option argument %s\n", m->argv[i]);
  }

  // Check if grid size is divisible by block size
  if ((grid_x % block_x) || (grid_y % block_y) || (grid_z % block_z)) {
    CkAbort("Grid size is not divisible by block size");
  }

  // Check if load balancing frequency is a multiple of refine frequency
  if (lb_freq % refine_freq) {
    CkAbort("LB frequency (%d) should be a multiple of refine requency (%d)",
        lb_freq, refine_freq);
  }

  // Set number of chares per dimension
  n_chares_x = grid_x / block_x;
  n_chares_y = grid_y / block_y;
  n_chares_z = grid_z / block_z;
  int n_chares = n_chares_x * n_chares_y * n_chares_z;

  // Set minimum depth
  // Only round it up if it is very close to the rounded up integer, otherwise
  // round it down
  float fdepth = log(n_chares) / log(NUM_CHILDREN);
  min_depth = (fabs(fdepth - ceil(fdepth)) < 0.000001) ? ceil(fdepth) : floor(fdepth);
  if (min_depth > max_depth) {
    CkAbort("Minimum depth > maximum depth: try increasing maximum depth");
  }
  CkAssert(min_depth > 0);

  // Initialize constants
  x_min = 0; x_max = 1;
  y_min = 0; y_max = 1;
  z_min = 0; z_max = 1;

  dx = (x_max - x_min) / float(grid_x);
  dy = (y_max - y_min) / float(grid_y);
  dz = (z_max - z_min) / float(grid_z);

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
           "* Grid dimensions: %d x %d x %d\n"
           "* Block dimensions: %d x %d x %d\n"
           "* Initial chares: %d x %d x %d\n"
           "* Minimum depth: %d\n"
           "* Maximum depth: %d\n"
           "* Number of iterations: %d\n"
           "* Refinement frequency: %d\n"
           "* Load balancing frequency: %d\n\n",
           grid_x, grid_y, grid_z,
           block_x, block_y, block_z,
           n_chares_x, n_chares_y, n_chares_z,
           min_depth, max_depth, max_iters, refine_freq, lb_freq);

  // Create tree of chares
  CProxy_AdvMap map = CProxy_AdvMap::ckNew();
  CkArrayOptions opts;
  opts.setMap(map);
  mesh = CProxy_MeshBlock::ckNew(opts);
  for (int i = 0; i < n_chares_x; ++i) {
    for (int j = 0; j < n_chares_y; ++j) {
      for (int k = 0; k < n_chares_z; ++k) {
        mesh[OctIndex(i, j, k, min_depth)].insert(x_min, x_max, y_min, y_max, z_min, z_max);
      }
    }
  }
  mesh.doneInserting();

  // Begin simulation
  CkStartQD(CkCallback(CkIndex_Main::startMeshGeneration(), thisProxy));
  mesh_manager = CProxy_MeshManager::ckNew();
}

void Main::startMeshGeneration() {
  start_time = CkWallTimer();
  mesh.iterate();
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
