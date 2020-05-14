#include "Headers.h"
#include <assert.h>

#define inInitialMeshGenerationPhase (gen_iter <= max_depth)

extern CProxy_Main main_proxy;
extern CProxy_MeshBlock mesh;
extern int grid_y, grid_x, grid_z;
extern int block_x, block_y, block_z;
extern int n_chares_x, n_chares_y, n_chares_z;
extern int min_depth, max_depth;
extern int max_iters, refine_freq, lb_freq;
extern bool verbose;
extern float vx, vy, vz;
extern float apx, anx, apy, any, apz, anz;
extern float x_ctr, y_ctr, z_ctr, radius;
extern float tmax, t, dt, cfl;
float refine_filter = 0.01;
float refine_cutoff= 0.2, derefine_cutoff = 0.05;

CProxy_MeshManager mesh_manager;

#ifdef USE_GPU
extern void gpuHostAlloc(void** ptr, size_t size);
extern void gpuHostFree(void* ptr);
extern void gpuDeviceAlloc(void** ptr, size_t size);
extern void gpuDeviceFree(void* ptr);
extern void gpuStreamCreate(cudaStream_t* stream_ptr);
extern void gpuStreamDestroy(cudaStream_t stream);

extern float invokeDecisionKernel(cudaStream_t, float*, float*, float*, float*,
    float*, float*, float, float, float, float, int, void*);
extern void invokeComputeKernel(cudaStream_t, float*, float*, float*, float*,
    float, float, float, float, float, float, float, float, float, float, int, void*);
#endif

enum BitsToUse { LOW, HIGH, BOTH };

static void populateIndices(BitsToUse x, BitsToUse y, BitsToUse z, std::vector<int>& ids)
{
  BitsToUse dims[3] = {z, y, x};

  for (unsigned i = 0; i < NUM_CHILDREN; ++i) {
    bool keep = true;
    for (unsigned j = 0; j < numDims; ++j) {
      unsigned bit = 1 << j;
      bool set = (i & bit) != 0;

      if (dims[j] == BOTH) continue;

      if (set == (dims[j] == LOW)) {
        keep = false;
        break;
      }
    }
    if (keep)
      ids.push_back(i);
  }
}

// Given an index, return the children that share the surface in the given direction
inline static void getChildrenInDir(OctIndex my_index, int dir, std::vector<OctIndex>& children)
{
  std::vector<int> ids;

  switch (dir) {
    case LEFT:     populateIndices(LOW,  BOTH, BOTH, ids); break;
    case RIGHT:    populateIndices(HIGH, BOTH, BOTH, ids); break;
    case UP:       populateIndices(BOTH, HIGH, BOTH, ids); break;
    case DOWN:     populateIndices(BOTH, LOW,  BOTH, ids); break;
    case FORWARD:  populateIndices(BOTH, BOTH, HIGH, ids); break;
    case BACKWARD: populateIndices(BOTH, BOTH, LOW,  ids); break;
  }

  for (auto iter = ids.begin(), end = ids.end(); iter != end; ++iter) {
    children.push_back(my_index.getChild(*iter));
  }
  CkAssert(children.size() == 4);
}

inline static void setFirstHalf(int& min, int& max) { max /= 2; }
inline static void setSecondHalf(int& min, int& max) { min = (max / 2) + 1; }

inline static void populateQuadrant(bool bit1, bool bit2, int& min1, int& max1, int& min2, int& max2)
{
  if (bit1) setSecondHalf(min1, max1);
  else      setFirstHalf(min1, max1);

  if (bit2) setSecondHalf(min2, max2);
  else      setFirstHalf(min2, max2);
}

enum {
  X_MASK = 1 << 2,
  Y_MASK = 1 << 1,
  Z_MASK = 1 << 0
};

inline static void populateRanges(int dir, int octant, int &x_min, int &x_max,
                           int& y_min, int& y_max, int& z_min, int& z_max)
{
  x_min = 0; x_max = block_x-1;
  y_min = 0; y_max = block_y-1;
  z_min = 0; z_max = block_z-1;

  switch (dir) {
    case UP:
      y_min = y_max = block_y+1;
      if (octant >= 0) {
        populateQuadrant(octant & X_MASK, octant & Z_MASK, x_min, x_max, z_min, z_max);
      }
      break;
    case DOWN:
      y_max = 0;
      if (octant >= 0) {
        populateQuadrant(octant & X_MASK, octant & Z_MASK, x_min, x_max, z_min, z_max);
      }
      break;
    case LEFT:
      x_max = 0;
      if (octant >= 0) {
        populateQuadrant(octant & Y_MASK, octant & Z_MASK, y_min, y_max, z_min, z_max);
      }
      break;
    case RIGHT:
      x_min = x_max = block_x+1;
      if (octant >= 0) {
        populateQuadrant(octant & Y_MASK, octant & Z_MASK, y_min, y_max, z_min, z_max);
      }
      break;
    case FORWARD:
      z_min = z_max = block_z+1;
      if (octant >= 0) {
        populateQuadrant(octant & X_MASK, octant & Y_MASK, x_min, x_max, y_min, y_max);
      }
      break;
    case BACKWARD:
      z_max = 0;
      if (octant >= 0) {
        populateQuadrant(octant & X_MASK, octant & Y_MASK, x_min, x_max, y_min, y_max);
      }
      break;
  }
}

inline static bool getOctantRange(int octant, int& x_min, int& x_max, int& y_min, int& y_max,
                           int& z_min, int& z_max)
{
  x_min = 1; x_max = block_x;
  y_min = 1; y_max = block_y;
  z_min = 1; z_max = block_z;

  if (octant & X_MASK)  x_min = block_x/2+1;
  else                  x_max = block_x/2;
  if (octant & Y_MASK)  y_min = block_y/2+1;
  else                  y_max = block_y/2;
  if (octant & Z_MASK)  z_min = block_z/2+1;
  else                  z_max = block_z/2;
}

MeshManager::MeshManager() : workUnitCount(0), compute_time_sum(0.0),
  compute_time_cnt(0), decision_time_sum(0.0), decision_time_cnt(0)
{
  // delu and delua are 4D arrays
#ifdef USE_GPU
  size_t delu_size = sizeof(float)*numDims*(block_x+2)*(block_y+2)*(block_z+2);
  gpuDeviceAlloc((void**)&d_delu, delu_size);
  gpuDeviceAlloc((void**)&d_delua, delu_size);
#endif
  delu = new float***[numDims];
  delua = new float***[numDims];

  for (int d = 0; d < numDims; ++d) {
    delu[d] = new float**[block_x+2];
    delua[d] = new float**[block_x+2];
    for (int i = 0; i < block_x+2; ++i) {
      delu[d][i] = new float*[block_y+2];
      delua[d][i] = new float*[block_y+2];
      for (int j = 0; j < block_y+2; ++j) {
        delu[d][i][j] = new float[block_z+2];
        delua[d][i][j] = new float[block_z+2];
      }
    }
  }
}

MeshManager::MeshManager(CkMigrateMessage* m) : CBase_MeshManager(m)
{
#ifdef USE_GPU
  size_t delu_size = sizeof(float)*numDims*(block_x+2)*(block_y+2)*(block_z+2);
  gpuDeviceAlloc((void**)&d_delu, delu_size);
  gpuDeviceAlloc((void**)&d_delua, delu_size);
#endif
  delu = new float***[numDims];
  delua = new float***[numDims];

  for (int d = 0; d < numDims; ++d) {
    delu[d] = new float**[block_x+2];
    delua[d] = new float**[block_x+2];
    for (int i = 0; i < block_x+2; ++i) {
      delu[d][i] = new float*[block_y+2];
      delua[d][i] = new float*[block_y+2];
      for (int j = 0; j < block_y+2; ++j) {
        delu[d][i][j] = new float[block_z+2];
        delua[d][i][j] = new float[block_z+2];
      }
    }
  }
}

MeshManager::~MeshManager()
{
#ifdef USE_GPU
  gpuDeviceFree(d_delu);
  gpuDeviceFree(d_delua);
#endif
}

void MeshManager::incrementWorkUnitCount(int iter)
{
  workUnitCount++;
  workUnits[iter]++;

  if (iter % lb_freq == 0 || iter % lb_freq == lb_freq-1) {
    // This is either a load balancing iteration or the one right before it
    minLoad[iter] += 1;
    maxLoad[iter] += 1;
    avgLoad[iter] += 1;
  }
}

void MeshManager::reduceWorkUnits()
{
  CkCallback cb(CkReductionTarget(Main,totalWorkUnits), main_proxy);
  contribute(sizeof(int), &workUnitCount, CkReduction::sum_int, cb);
}

void MeshManager::processQdTimes(std::map<int, pair<float, float>> peQdtimes,
                                    std::map<int, pair<float, float>> peRemeshtimes,
                                    std::map<int, int> peWorkunits, std::map<int, int> peminLoad,
                                    std::map<int, int> pemaxLoad, std::map<int, float> peavgLoad)
{
  for (auto it = peQdtimes.begin(); it != peQdtimes.end(); it++) {
    if (qdtimes.find(it->first) == qdtimes.end()) {
      qdtimes[it->first] = std::pair<float, float>(std::numeric_limits<float>::min(),
                                                   std::numeric_limits<float>::max());
    }
    qdtimes[it->first].first = std::max(qdtimes[it->first].first, it->second.first);
    qdtimes[it->first].second = std::min(qdtimes[it->first].second, it->second.second);
  }

  for (auto it = peRemeshtimes.begin(); it != peRemeshtimes.end(); it++) {
    if (remeshtimes.find(it->first) == remeshtimes.end()) {
      remeshtimes[it->first] = std::pair<float, float>(0, std::numeric_limits<float>::max());
    }
    remeshtimes[it->first].first = std::max(remeshtimes[it->first].first, it->second.first);
    remeshtimes[it->first].second = std::min(remeshtimes[it->first].second, it->second.second);
  }

  for (auto it = peWorkunits.begin(); it != peWorkunits.end(); it++) {
    workUnits[it->first] += it->second;
    if (it->first % lb_freq == 0 || it->first % lb_freq == lb_freq-1) {
      // This is either a load balancing iteration or the one right before it
      minLoad[it->first] = std::min(minLoad[it->first], peminLoad[it->first]);
      maxLoad[it->first] = std::max(maxLoad[it->first], pemaxLoad[it->first]);
      avgLoad[it->first] += peavgLoad[it->first];
    }
  }
}

void MeshManager::printLogs()
{
  ckout << "Compute function average time: " << compute_time_sum/compute_time_cnt << endl;
  ckout << "Refinement decision average time: " << decision_time_sum/decision_time_cnt << endl;
  ckout << endl;

  if (verbose) {
    ckout << "QD times (iteration, duration):" << endl;
    for (auto it = qdtimes.begin(); it != qdtimes.end(); it++) {
      ckout << "(" << it->first << ", " << it->second.second - it->second.first << ")\n";
    }
    ckout << endl;

    ckout << "Remeshing times (iteration, duration):" << endl;
    for (auto it = remeshtimes.begin(); it != remeshtimes.end(); it++) {
      ckout << "(" << it->first << ", " << it->second.second - it->second.first << ")\n";
    }
    ckout << endl;

    ckout << "Work units per iteration (iteration, count):" << endl;
    for (auto it = workUnits.begin(); it != workUnits.end(); it++) {
      ckout << "(" << it->first << ", " << it->second << ")\n";
    }
    ckout << endl;

    ckout << "Load balancing stats (iteration, average load, min load, max load):" << endl;
    for (auto it = minLoad.begin(); it != minLoad.end(); it++) {
      avgLoad[it->first] /= CkNumPes();
      ckout << "(" << it->first << ", " << int(avgLoad[it->first]*100)/100.
        << ", " << minLoad[it->first] << ", " << maxLoad[it->first] << ")\n";
    }
    ckout << endl;
  }

  CkExit();
}

void MeshManager::meshGenerationPhaseIsOver() {}

void MeshBlock::packGhosts(){
  recv_count = 0;
  recv_ghosts.clear();

  for (auto it = neighbors.begin(), iend = neighbors.end(); it != iend; ++it) {
    it->second.setDataSent(false);
  }

  // Pack ghosts into contiguous blocks
  // YZ surfaces
  for (int j = 1; j <= block_y; ++j) {
    for (int k = 1; k <= block_z; ++k) {
      left_surface[index_yz(j-1,k-1)] = u[index(1,j,k)];
      right_surface[index_yz(j-1,k-1)] = u[index(block_x,j,k)];
    }
  }
  // XZ Surfaces
  for (int i = 1; i <= block_x; ++i) {
    for (int k = 1; k <= block_z; ++k) {
      top_surface[index_xz(i-1,k-1)] = u[index(i, block_y, k)];
      bottom_surface[index_xz(i-1,k-1)] = u[index(i, 1, k)];
    }
  }
  // XY surfaces
  for (int i = 1; i <= block_x; ++i) {
    for (int j = 1; j <= block_y; ++j) {
      forward_surface[index_xy(i-1,j-1)] = u[index(i, j, block_z)];
      backward_surface[index_xy(i-1,j-1)] = u[index(i, j, 1)];
    }
  }
}

void MeshBlock::applyInitialCondition(){
  VB(logfile << "applying initial condition" << std::endl;);
  float rcub;
  for(int i=0; i<block_x+2; i++)
    for(int j=0; j<block_y+2; j++)
      for(int k=0; k<block_z+2; k++){
        rcub = (x[i] - x_ctr)*(x[i]-x_ctr) +
               (y[j] - y_ctr)*(y[j]-y_ctr) +
               (z[k] - z_ctr)*(z[k]-z_ctr);
        u[index(i, j, k)] = (rcub<=radius*radius) ? 2:1;
        //VB(logfile << x[i] << " " << y[j] << " " << z[k] << ": " << rcub << " " << radius*radius << " " << u[index(i, j, k)] << std::endl;);
      }
}

inline void MeshBlock::hostAlloc(float* &p, int size){
  p = new float[size];
}

inline void MeshBlock::hostFree(float* p) {
  delete [] p;
}

void MeshBlock::allAlloc(){
#ifdef USE_GPU
  gpuHostAlloc((void**)&u, (block_x+2)*(block_y+2)*(block_z+2)*sizeof(float));
  gpuHostAlloc((void**)&u2, (block_x+2)*(block_y+2)*(block_z+2)*sizeof(float));
  gpuHostAlloc((void**)&u3, (block_x+2)*(block_y+2)*(block_z+2)*sizeof(float));
  gpuHostAlloc((void**)&h_error, sizeof(float));
  gpuDeviceAlloc((void**)&d_u, sizeof(float)*(block_x+2)*(block_y+2)*(block_z+2));
  gpuDeviceAlloc((void**)&d_u2, sizeof(float)*(block_x+2)*(block_y+2)*(block_z+2));
  gpuDeviceAlloc((void**)&d_u3, sizeof(float)*(block_x+2)*(block_y+2)*(block_z+2));
  gpuDeviceAlloc((void**)&d_error, sizeof(float));

  gpuStreamCreate(&computeStream);
  gpuStreamCreate(&decisionStream);
#else
  hostAlloc(u, (block_x+2)*(block_y+2)*(block_z+2));
  hostAlloc(u2, (block_x+2)*(block_y+2)*(block_z+2));
  hostAlloc(u3, (block_x+2)*(block_y+2)*(block_z+2));
#endif

  hostAlloc(x, block_x+2);
  hostAlloc(y, block_y+2);
  hostAlloc(z, block_z+2);

  hostAlloc(left_surface, block_y*block_z);
  hostAlloc(right_surface, block_y*block_z);
  hostAlloc(top_surface, block_x*block_z);
  hostAlloc(bottom_surface, block_x*block_z);
  hostAlloc(forward_surface, block_x*block_y);
  hostAlloc(backward_surface, block_x*block_y);
}

void MeshBlock::allFree(){
#ifdef USE_GPU
  gpuHostFree(u);
  gpuHostFree(u2);
  gpuHostFree(u3);
  gpuHostFree(h_error);
  gpuDeviceFree(d_u);
  gpuDeviceFree(d_u2);
  gpuDeviceFree(d_u3);
  gpuDeviceFree(d_error);

  gpuStreamDestroy(computeStream);
  gpuStreamDestroy(decisionStream);
#else
  hostFree(u);
  hostFree(u2);
  hostFree(u3);
#endif

  hostFree(x);
  hostFree(y);
  hostFree(z);

  hostFree(left_surface);
  hostFree(right_surface);
  hostFree(top_surface);
  hostFree(bottom_surface);
  hostFree(forward_surface);
  hostFree(backward_surface);
}

void MeshBlock::init() {
  mesh_manager_local = mesh_manager.ckLocalBranch();

  allAlloc();

  usesAutoMeasure = false;
  usesAtSync = true;
  remeshStartTime = 0;

  if(inInitialMeshGenerationPhase){
    for(int i=0; i<block_x+2; i++)
      x[i] = x_min + float(i)*dx - 0.5*dx;


    for(int i=0; i<block_y+2; i++)
      y[i] = y_min + float(i)*dy - 0.5*dy;

    for(int i=0; i<block_z+2; i++)
      z[i] = z_min + float(i)*dz - 0.5*dz;

    applyInitialCondition();
  }

  FOR_EACH_CHILD
    child_isRefined[i] = false;
  END_FOR
  this->isRefined = false;

  FOR_EACH_NEIGHBOR
    neighbors[thisIndex.getNeighbor(i)] = Neighbor(i);
  END_FOR

  parent = (thisIndex.getDepth()>min_depth) ? thisIndex.getParent() : thisIndex;
  resetMeshRestructureData();
  phase1Over = false;
}

// Constructor used to create initial chares
MeshBlock::MeshBlock(float x_min, float x_max, float y_min, float y_max,
                     float z_min, float z_max)
{
  thisIndex.getCoordinates(xc, yc, zc);

  dx = (x_max - x_min) / float(grid_x);
  dy = (y_max - y_min) / float(grid_y);
  dz = (z_max - z_min) / float(grid_z);

  nx = grid_x / n_chares_x;
  ny = grid_y / n_chares_y;
  nz = grid_z / n_chares_z;

  myt = t;
  mydt = dt;

  this->x_min = xc * nx * dx;
  this->y_min = yc * ny * dy;
  this->z_min = zc * nz * dz;

  iter = 0;
  gen_iter = 0;

  init();
}

// Constructor used to create children chares on refinement
MeshBlock::MeshBlock(float dx, float dy, float dz, float myt, float mydt,
                     float x_min, float y_min, float z_min,
                     int gen_iter, int iter, std::vector<float> refined_u,
                     std::map<OctIndex, Neighbor> parent_neighbors)
{
  thisIndex.getCoordinates(xc, yc, zc);

  this->dx = dx;
  this->dy = dy;
  this->dz = dz;

  this->myt = myt;
  this->mydt = mydt;

  this->x_min = x_min;
  this->y_min = y_min;
  this->z_min = z_min;

  nx = grid_x / n_chares_x;
  ny = grid_y / n_chares_y;
  nz = grid_z / n_chares_z;

  this->gen_iter = gen_iter;
  this->iter = iter;

  init();

  auto it = neighbors.begin();
  while (it != neighbors.end()) {
    const OctIndex& neighbor_idx = it->first;
    Neighbor& neighbor = it->second;
    bool should_delete = false;

    if (neighbor_idx.getParent() != thisIndex.getParent()) {
      auto parent_it = parent_neighbors.find(neighbor_idx.getParent());
      if (parent_it == parent_neighbors.end()) {
        should_delete = true;
      } else {
        Neighbor& parent_neighbor = parent_it->second;
        if (parent_neighbor.isRefined()) {
          switch (parent_neighbor.getDecision(neighbor_idx.getOctant())) {
            case COARSEN: should_delete = true; break;
            case REFINE: neighbor.setRefined(true); break;
          }
        } else {
          switch (parent_neighbor.getDecision()) {
            case STAY: should_delete = true; break;
            case COARSEN: CkAbort("This neighbor cannot derefine");
          }
        }
      }
    }

    // Advance iterator, removing this index/neighbor pair if requested
    if (should_delete) {
      neighbors.erase(it++);
    } else {
      it++;
    }
  }
  CkAssert(neighbors.size() >= numDims);

  // Ensure 3 of the neighbors share the same parent as me
  unsigned same_parent = 0;
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    if (it->first.getParent() == thisIndex.getParent()) {
      same_parent++;
    }
  }
  CkAssert(same_parent == 3);

  // Fill in the refined data unless we are in the initial mesh generation phase
  if (!inInitialMeshGenerationPhase) {
    int elem = 0;
    for (int k = 1; k <= block_z; k++) {
      for (int j = 1; j <= block_y; j++) {
        for (int i = 1; i <= block_x; i++) {
          u[index(i,j,k)] = refined_u[elem++];
        }
      }
    }
  }

  // Begin iterating
  thisProxy[thisIndex].iterate();
}

void MeshBlock::pup(PUP::er &p) {
  p|isRefined;
  p|depth;

  PUParray(p, child_isRefined, NUM_CHILDREN);

  p|recv_ghosts;

  p|decision;
  PUParray(p, child_decision, NUM_CHILDREN);

  p|parentHasAlreadyMadeDecision;
  p|hasReceivedParentDecision;
  PUParray(p, fname, 100);
  p|parent;
  p|neighbors;
  p|uncleDecisions;

  p|xc;
  p|yc;
  p|zc;
  p|recv_count;
  p|phase1Over;

  if(p.isUnpacking()){
    allAlloc(); // Needed as there is no migration constructor
    //resetMeshRestructureData();
  }

  PUParray(p, u, (block_x+2)*(block_y+2)*(block_z+2));
  PUParray(p, x, block_x+2);
  PUParray(p, y, block_y+2);
  PUParray(p, z, block_z+2);

  p|iter;
  p|gen_iter;
  p|up_x;
  p|un_x;
  p|up_y;
  p|un_y;
  p|up_z;
  p|un_z;
  p|myt;
  p|mydt;

  p|dx; p|dy; p|dz;
  p|nx; p|ny; p|nz;
  p|x_min; p|x_max;
  p|y_min; p|y_max;
  p|z_min; p|z_max;

  p|itBeginTime;
  p|remeshStartTime;
  p|remeshEndTime;
  p|lastBusyTime; p|lastIdleTime;
  //p|logfile;
  p|lastIdleTimeQD;
  p|finishedPhase1;
  p|phase1Over;
  p|amr3d_i;
  p|ichild;
}

MeshBlock::~MeshBlock(){
  if (isLeaf) {
    allFree(); // FIXME: Shouldn't this be called for ALL chares?
  }
}

inline float downSample(float* u, int x, int y, int z) {
  return (
    u[index(x, y, z)] + u[index(x+1, y, z)] +
    u[index(x, y+1, z)] + u[index(x+1, y+1, z)] +
    u[index(x, y, z+1)] + u[index(x+1, y, z+1)] +
    u[index(x, y+1, z+1)] + u[index(x+1, y+1, z+1)]) / 8.0;
}

class SurfaceIterator {
  float *u;
  int minx, maxx;
  int miny, maxy;
  int minz, maxz;
  int curx, cury, curz;
  int dx, dy, dz;
  bool done;

  void incX() {
    curx += dx;
    if (curx > maxx || (dx == 0)) {
      curx = minx;
      incY();
    }
  }

  void incY() {
    cury += dy;
    if (cury > maxy || (dy == 0)) {
      cury = miny;
      incZ();
    }
  }

  void incZ() {
    curz += dz;
    if (curz > maxz || (dz == 0)) {
      done = true;
    }
  }

public:
  SurfaceIterator(float *u,
                   int minx, int maxx, int dx,
                   int miny, int maxy, int dy,
                   int minz, int maxz, int dz)
  : u(u), minx(minx), maxx(maxx), dx(dx),
          miny(miny), maxy(maxy), dy(dy),
          minz(minz), maxz(maxz), dz(dz),
          curx(minx), cury(miny), curz(minz),
          done(false) {
    assert((dx == 0) || (minx < maxx));
    assert((dy == 0) || (miny < maxy));
    assert((dz == 0) || (minz < maxz));
  }

  int getX() { assert(!done); return curx; }
  int getY() { assert(!done); return cury; }
  int getZ() { assert(!done); return curz; }

  inline SurfaceIterator& operator++() {
    assert(!done);
    incX();
    return *this;
  }

  inline float left()     { assert(!done); return u[index(curx-1, cury, curz)];  }
  inline float right()    { assert(!done); return u[index(curx+1, cury, curz)];  }
  inline float up()       { assert(!done); return u[index(curx, cury+1, curz)];  }
  inline float down()     { assert(!done); return u[index(curx, cury-1, curz)];  }
  inline float forward()  { assert(!done); return u[index(curx, cury, curz+1)];  }
  inline float backward() { assert(!done); return u[index(curx, cury, curz-1)];  }

  inline float& operator*() {
    assert(!done);
    return u[index(curx, cury, curz)];
  }

  bool isDone() {
    return done;
  }
};

inline float* MeshBlock::getGhostBuffer(int dir) {
  switch (dir) {
    case UP:       return top_surface;
    case DOWN:     return bottom_surface;
    case LEFT:     return left_surface;
    case RIGHT:    return right_surface;
    case FORWARD:  return forward_surface;
    case BACKWARD: return backward_surface;
  }
}

inline int MeshBlock::getGhostSize(int dir) {
  switch (dir) {
    case UP:
    case DOWN:
      return block_x * block_z;
    case RIGHT:
    case LEFT:
      return block_y * block_z;
    case FORWARD:
    case BACKWARD:
      return block_x * block_y;
  }
}

inline int MeshBlock::getSourceDirection(int dir) {
  switch (dir) {
    case UP:       return DOWN;
    case DOWN:     return UP;
    case LEFT:     return RIGHT;
    case RIGHT:    return LEFT;
    case FORWARD:  return BACKWARD;
    case BACKWARD: return FORWARD;
  }
}

void MeshBlock::sendGhost(int dir){
  int ghost_size = getGhostSize(dir);
  float* ghost_data = NULL;

  OctIndex neighbor_idx = thisIndex.getNeighbor(dir);
  auto it = neighbors.find(neighbor_idx);

  if (it == neighbors.end()) {
    // Uncle case, neighbor doesn't exist in this direction (at this level)
    OctIndex receiver = neighbor_idx.getParent();

    // TODO: Remove use of subdirections
    ghost_data = getGhostBuffer(dir);
    ghost_size /= 4;

    int x_min, x_max, y_min, y_max, z_min, z_max;
    populateRanges(dir, -1, x_min, x_max, y_min, y_max, z_min, z_max);

    int dx = (x_min == x_max) ? 0 : 2;
    int dy = (y_min == y_max) ? 0 : 2;
    int dz = (z_min == z_max) ? 0 : 2;

    if (x_min != x_max) x_min++;
    else if (x_min == 0) x_min = x_max = 1;
    else if (x_min == block_x+1) x_min = x_max = block_x-1;

    if (y_min != y_max) y_min++;
    else if (y_min == 0) y_min = y_max = 1;
    else if (y_min == block_y+1) y_min = y_max = block_y-1;

    if (z_min != z_max) z_min++;
    else if (z_min == 0) z_min = z_max = 1;
    else if (z_min == block_z+1) z_min = z_max = block_z-1;

    // Iterate the original surface and downsample to fill in the ghost surface
    SurfaceIterator surface_it(u, x_min, x_max, dx,
                               y_min, y_max, dy,
                               z_min, z_max, dz);
    int elem;
    for (elem = 0; !surface_it.isDone(); ++it, ++elem) {
      ghost_data[elem] = downSample(u, surface_it.getX(), surface_it.getY(),
                                    surface_it.getZ());
    }
    CkAssert(elem == ghost_size);

    // Send the downsampled ghost data
    thisProxy(receiver).receiveGhost(iter, getSourceDirection(dir),
                                      thisIndex.getOctant(), ghost_size,
                                      ghost_data);
  } else if (!it->second.isRefined()) {
    // Friend at same level
    ghost_data = getGhostBuffer(dir);
    thisProxy(neighbor_idx).receiveGhost(iter, getSourceDirection(dir), -1,
                                          ghost_size, ghost_data);
  }
}

void MeshBlock::processGhost(int iter, int dir, int quadrant, int size, float gh[]) {
  bool from_nephew = (quadrant >= 0);

  OctIndex neighbor_idx = thisIndex.getNeighbor(dir);
  if (from_nephew) neighbor_idx = neighbor_idx.getChild(quadrant);
  recv_ghosts.insert(neighbor_idx);
  recv_count += (from_nephew) ? 0.25 : 1;

  int x_min, x_max, y_min, y_max, z_min, z_max;
  populateRanges(dir, quadrant, x_min, x_max, y_min, y_max, z_min, z_max);

  int dx = (x_min == x_max) ? 0 : 1;
  int dy = (y_min == y_max) ? 0 : 1;
  int dz = (z_min == z_max) ? 0 : 1;

  if (x_max != x_min) {
    x_min++;
    x_max++;
  }
  if (y_max != y_min) {
    y_min++;
    y_max++;
  }
  if (z_max != z_min) {
    z_min++;
    z_max++;
  }

  // Iterate the received ghost surface and fill in the original surface
  SurfaceIterator surface_it(u, x_min, x_max, dx,
                             y_min, y_max, dy,
                             z_min, z_max, dz);
  for (int i = 0; i < size; ++i, ++surface_it) {
    *surface_it = gh[i];
  }
  CkAssert(surface_it.isDone());
}

void MeshBlock::sendReadyData(){
  // Check if data can be sent to any of the refined neighbors.
  // If the neighbors are at the same level or do not exist at all,
  // data will be sent in begin_iteration and need not be sent here.
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    Neighbor& neighbor = it->second;
    if (neighbor.isRefined() && !neighbor.isDataSent()) {
      // Find neighboring children
      std::vector<OctIndex> children;
      getChildrenInDir(thisIndex, neighbor.getDir(), children);

      std::vector<OctIndex> all_neighbors;
      for (auto child_it = children.begin(); child_it != children.end(); ++child_it) {
        std::vector<OctIndex> child_neighbors = child_it->getNeighbors();
        all_neighbors.insert(all_neighbors.end(),
                             child_neighbors.begin(), child_neighbors.end());
      }

      bool received_all = true;
      for (auto ait = all_neighbors.begin(); ait != all_neighbors.end(); ++ait) {
        OctIndex& neighbor_idx = *ait;
        OctIndex neighbor_parent = neighbor_idx.getParent();
        if (neighbor_parent != thisIndex) { // External neighbor
          if (!recv_ghosts.count(neighbor_idx) && !recv_ghosts.count(neighbor_parent)) {
            received_all = false;
            break;
          }
        }
      }

      if (received_all) {
        interpolateAndSend(neighbor.getDir());
        neighbor.setDataSent(true);
      }
    }
  }
}

void MeshBlock::interpolateAndSend(int dir) {
  int uncledir = getSourceDirection(dir);
  std::vector<OctIndex> children;
  getChildrenInDir(thisIndex.getNeighbor(dir), uncledir, children);
  for (std::vector<OctIndex>::iterator I = children.begin(),
       E = children.end(); I != E; ++I)
    interpolateAndSendToNephew(uncledir, *I);
}

void MeshBlock::interpolateAndSendToNephew(int uncledir, OctIndex QI) {
  float *boundary;
  int count = getGhostSize(uncledir);

  // TODO: Better generalize this code!
  int p = 1, m = -1;
  float sx_r, sx_l, sy_u, sy_d, sz_f, sz_b;

  float *sx[]={&sx_l, &sx_r};
  float *sy[]={&sy_d, &sy_u};
  float *sz[]={&sz_b, &sz_f};

  int octant = QI.getOctant();

  int x_min, x_max, y_min, y_max, z_min, z_max;
  populateRanges(getSourceDirection(uncledir), octant,
                 x_min, x_max,
                 y_min, y_max,
                 z_min, z_max);

  int dx = 1, dy = 1, dz = 1;
  unsigned columncount;
  float **a, **b, *c;
  if (x_min == x_max) {
    a = sy;
    b = sz;
    c = (x_max == 0) ? &sx_l : &sx_r;
    columncount = block_y;
    dx = 0;
    y_min++; y_max++;
    z_min++; z_max++;
    if(x_min==0){
        x_min = x_max = 1;
    }else{
        x_min = x_max = block_x;
    }

  }
  if (y_min == y_max) {
    a = sx;
    b = sz;
    c = (y_max == 0) ? &sy_d : &sy_u;
    columncount = block_x;
    dy = 0;
    x_min++; x_max++;
    z_min++; z_max++;
    if(y_min==0){
        y_min = y_max = 1;
    }else{
        y_min = y_max = block_y;
    }
  }
  if (z_min == z_max) {
    a = sx;
    b = sy;
    c = (z_max == 0) ? &sz_b : &sz_f;
    columncount = block_x;
    dz = 0;
    x_min++; x_max++;
    y_min++; y_max++;
    if(z_min==0){
        z_min = z_max = 1;
    }else{
        z_min = z_max = block_z;
    }
  }

  SurfaceIterator in(u, x_min, x_max, dx,
                         y_min, y_max, dy,
                         z_min, z_max, dz);

  boundary = getGhostBuffer(getSourceDirection(uncledir));

  unsigned counter = 0;
  for (; !in.isDone(); ++in) {
    sx_r = (in.right()      - *in) / 4;
    sx_l = (in.left()       - *in) / 4;
    sy_u = (in.up()         - *in) / 4;
    sy_d = (in.down()       - *in) / 4;
    sz_f = (in.forward()    - *in) / 4;
    sz_b = (in.backward()  -  *in) / 4;

    int a_pos = counter % columncount;
    int b_pos = counter / columncount;
    counter++;

    boundary[index2d(2*a_pos,   2*b_pos,   columncount)] = *in + *a[0] + *b[0] + *c;
    boundary[index2d(2*a_pos+1, 2*b_pos,   columncount)] = *in + *a[1] + *b[0] + *c;
    boundary[index2d(2*a_pos,   2*b_pos+1, columncount)] = *in + *a[0] + *b[1] + *c;
    boundary[index2d(2*a_pos+1, 2*b_pos+1, columncount)] = *in + *a[1] + *b[1] + *c;
  }
  assert(counter == count/(1 << (numDims-1)));

  thisProxy(QI).receiveGhost(iter, uncledir, -2, count, boundary);
}

void MeshBlock::computeDone() {
#if defined(USE_GPU) && defined(USE_HAPI)
  double time_dur = CkWallTimer() - compute_start_time;
  mesh_manager_local->addComputeTime(time_dur);

  iterate();
#endif
}

void MeshBlock::compute(){
  //if(iter==1){
#ifdef LOGGER
    char logfilename[100];
    sprintf(logfilename, "out/snap_%s_%d.vtk", thisIndex.getIndexString().c_str(), iter);
    logfile.open(logfilename);
    logfile.precision(8);
    //logfile << "Variables=\"X\",\"Y\",\"Z\",\"radius\",\"nodeid\"\n";
    //logfile << "Zone I = 16, J = 16, K = 16, F = POINT\n";
    printData();
    logfile.close();
    /*int dims[] = {block_y+1, block_x+1, block_z+1};
    int vardims[] = {1};
    int centering[] = {0};
    const char *varnames[]  = {"viscosity"};
    //varnames = new char*[1];
    //varnames[0] = "viscosity";
    
    float *vars[1];
    vars[0] = new float[block_x*block_y*block_z];
    for(int i=0; i<block_x; i++)
        for(int j=0; j<block_y; j++)
            for(int k=0; k<block_z; k++){
                vars[0][(k*block_y+j)*block_x+i] = u[index(i+1,j+1,k+1)];
                //ckout << x[i+1] << " " << y[j+1] << " " << z[k+1] << " " << u[index(i+1,j+1,k+1)] << endl;
            }
    float *xn = new float[block_x+1];
    float *yn = new float[block_y+1];
    float *zn = new float[block_z+1];
    for(int i=0; i<block_x+1; i++) {xn[i]=x_min+dx*i;}
    for(int i=0; i<block_y+1; i++){yn[i]=y_min+dy*i;}
    for(int i=0; i<block_z+1; i++) {zn[i]=z_min+dz*i;}
    //float *vars[] = {u};
    write_rectilinear_mesh(logfilename, 0, dims, xn, yn, zn, 1, vardims, centering, varnames, vars);
    delete [] vars[0]; delete [] xn; delete [] yn; delete [] zn;
    for(int k=0; k<=block_z+1; k++){
        for(int j=0; j<=block_y+1; j++){
            for(int i=0; i<=block_x+1; i++){
                if(y[j]==0.484375)
                    logfile << std::setw(10) << u[index(i,j,k)] << " ";
            }
        }
        logfile << std::endl;
    }*/
#endif
  //}

  compute_start_time = CkWallTimer();

#ifndef USE_GPU
  /********** CPU CODE **********/
  memcpy(u2, u, sizeof(float)*(block_x+2)*(block_y+2)*(block_z+2));
  memcpy(u3, u, sizeof(float)*(block_x+2)*(block_y+2)*(block_z+2));
  FOR_EACH_ZONE
      up_x = (u[index(i+1,j,k)] - u[index(i,j,k)])/dx;
      un_x = (u[index(i,j,k)] - u[index(i-1,j,k)])/dx;
      up_y = (u[index(i,j+1,k)] - u[index(i,j,k)])/dy;
      un_y = (u[index(i,j,k)] - u[index(i,j-1,k)])/dy;
      up_z = (u[index(i,j,k+1)] - u[index(i,j,k)])/dz;
      un_z = (u[index(i,j,k)] - u[index(i,j,k-1)])/dz;

      u2[index(i,j,k)] = u[index(i,j,k)] - dt* (apx*un_x + anx*up_x) - dt*(apy*un_y + any*up_y) - dt*(apz*un_z + anz*up_z);
  END_FOR

  // =========================================
  FOR_EACH_ZONE
      up_x = (u2[index(i+1,j,k)] - u2[index(i,j,k)])/dx;
      un_x = (u2[index(i,j,k)] - u2[index(i-1,j,k)])/dx;
      up_y = (u2[index(i,j+1,k)] - u2[index(i,j,k)])/dy;
      un_y = (u2[index(i,j,k)] - u2[index(i,j-1,k)])/dy;
      up_z = (u2[index(i,j,k+1)] - u2[index(i,j,k)])/dz;
      un_z = (u2[index(i,j,k)] - u2[index(i,j,k-1)])/dz;

      u3[index(i,j,k)] = u2[index(i,j,k)] - dt* (apx*un_x + anx*up_x) - dt*(apy*un_y + any*up_y) - dt*(apz*un_z + anz*up_z);
      u[index(i,j,k)] = 0.5*(u2[index(i,j,k)] + u3[index(i,j,k)]);
  END_FOR
#else
  /********** GPU CODE **********/
#ifndef USE_HAPI
  invokeComputeKernel(computeStream, u, d_u, d_u2, d_u3, dx, dy, dz, dt, apx, apy, apz, anx, any, anz, block_x, NULL);
#else
  // create callback
  CkArrayIndexOctIndex myIndex = CkArrayIndexOctIndex(thisIndex);
  CkCallback *cb = new CkCallback(CkIndex_MeshBlock::computeDone(), myIndex, thisProxy);

  invokeComputeKernel(computeStream, u, d_u, d_u2, d_u3, dx, dy, dz, dt, apx, apy, apz, anx, any, anz, block_x, cb);
#endif // USE_HAPI
#endif // USE_GPU

#ifndef USE_HAPI
  double compute_time = CkWallTimer() - compute_start_time;
  mesh_manager_local->addComputeTime(compute_time);
#endif

#ifndef USE_HAPI
  iterate();
#endif
}

void MeshBlock::gotErrorFromGPU() {
#if defined(USE_GPU) && defined(USE_HAPI)
  double decision_time = CkWallTimer() - decision_start_time;
  mesh_manager_local->addDecisionTime(decision_time);

  float error = sqrt(*h_error);

  Decision newDecision;
  if (error < derefine_cutoff && thisIndex.getDepth() > min_depth) {
    newDecision = COARSEN;
  }
  else if (error > refine_cutoff && thisIndex.getDepth() < max_depth) {
    newDecision = REFINE;
  }
  else {
    newDecision = STAY;
  }

  newDecision = (decision != REFINE) ? std::max(decision, newDecision) : decision;
  VB(logfile << thisIndex.getIndexString().c_str() << " decision = " << newDecision << std::endl;);
  updateDecisionState(1, newDecision);
#endif
}

Decision MeshBlock::getGranularityDecision(){
  float delx = 0.5/dx;
  float dely = 0.5/dy;
  float delz = 0.5/dz;
  float error=0;

  decision_start_time = CkWallTimer();

//#ifndef USE_GPU
#if 1 // FIXME Do not use the GPU kernels for now, the error values are slightly different
  /********** CPU CODE **********/
  for(int i=1; i <= block_x; i++){
    for(int j=1; j<=block_y; j++){
      for(int k=1; k<=block_z; k++){
        // d/dx
        mesh_manager_local->delu[0][i][j][k] = (u[index(i+1, j, k)] - u[index(i-1, j, k)])*delx;
        mesh_manager_local->delua[0][i][j][k] = abs(u[index(i+1, j, k)]) + abs(u[index(i-1, j, k)])*delx;

        // d/dy
        mesh_manager_local->delu[1][i][j][k] = (u[index(i, j+1, k)] - u[index(i, j-1, k)])*dely;
        mesh_manager_local->delua[1][i][j][k] = (abs(u[index(i, j+1, k)]) + abs(u[index(i, j-1, k)]))*dely;

        // d/dz
        mesh_manager_local->delu[2][i][j][k] = (u[index(i, j, k+1)] - u[index(i, j, k-1)])*delz;
        mesh_manager_local->delua[2][i][j][k] = (abs(u[index(i, j, k+1)]) + abs(u[index(i, j, k-1)]))*delz;
      }
    }
  }

  int istart=2, iend=block_x-1,
      jstart=2, jend=block_y-1,
      kstart=2, kend=block_z-1;
  for (int i=istart;i<=iend;i++){
    for (int j=jstart;j<=jend;j++){
      for (int k=kstart;k<=kend;k++){
        for (int d = 0; d < numDims; ++d) {
          mesh_manager_local->delu2[3*d+0] = (mesh_manager_local->delu[d][i+1][j][k] - mesh_manager_local->delu[d][i-1][j][k])*delx;
          mesh_manager_local->delu3[3*d+0] = (abs(mesh_manager_local->delu[d][i+1][j][k]) + abs(mesh_manager_local->delu[d][i-1][j][k]))*delx;
          mesh_manager_local->delu4[3*d+0] = (mesh_manager_local->delua[d][i+1][j][k] + mesh_manager_local->delua[d][i-1][j][k])*delx;

          mesh_manager_local->delu2[3*d+1] = (mesh_manager_local->delu[d][i][j+1][k] - mesh_manager_local->delu[d][i][j-1][k])*dely;
          mesh_manager_local->delu3[3*d+1] = (abs(mesh_manager_local->delu[d][i][j+1][k]) + abs(mesh_manager_local->delu[d][i][j-1][k]))*dely;
          mesh_manager_local->delu4[3*d+1] = (mesh_manager_local->delua[d][i][j+1][k] + mesh_manager_local->delua[d][i][j-1][k])*dely;

          mesh_manager_local->delu2[3*d+2] = (mesh_manager_local->delu[d][i][j][k+1] - mesh_manager_local->delu[d][i][j][k-1])*delz;
          mesh_manager_local->delu3[3*d+2] = (abs(mesh_manager_local->delu[d][i][j][k+1]) + abs(mesh_manager_local->delu[d][i][j][k-1]))*delz;
          mesh_manager_local->delu4[3*d+2] = (mesh_manager_local->delua[d][i][j][k+1] + mesh_manager_local->delua[d][i][j][k-1])*delz;
        }

        // compute the error
        float num = 0.;
        float denom = 0.;

        for (int kk = 0; kk < numDims2; kk++){  // kk= 1, 2, 3, 4, 5, ... 9
          num = num + pow(mesh_manager_local->delu2[kk],2.);
          denom = denom + pow(mesh_manager_local->delu3[kk], 2.) + (refine_filter*mesh_manager_local->delu4[kk])*2;
        }
        // compare the square of the error
        if (denom == 0. && num != 0.){
          error = std::numeric_limits<float>::max();
        } else if (denom != 0.0){
          error = std::max(error, num/denom);
        }
      }
    }
  }

  double decision_time = CkWallTimer() - decision_start_time;
  mesh_manager_local->addDecisionTime(decision_time);

  error = sqrt(error);
  //CkPrintf("[Iter %d, Chare %d-%d-%d] error: %f\n", iter, xc, yc, zc, error);
  if(error < derefine_cutoff && thisIndex.getDepth() > min_depth) return COARSEN;
  else if(error > refine_cutoff && thisIndex.getDepth() < max_depth) return REFINE;
  else return STAY;
#else
  /********** GPU CODE **********/
//#ifndef USE_HAPI
#if 1 // FIXME Don't use HAPI version because it sometimes results in different refinement decisions
  // execute GPU kernel
  float error_gpu = invokeDecisionKernel(decisionStream, u, h_error, d_error, d_u, mesh_manager_local->d_delu, mesh_manager_local->d_delua, refine_filter, dx, dy, dz, block_x, NULL);

  double decision_time = CkWallTimer() - decision_start_time;
  mesh_manager_local->addDecisionTime(decision_time);

  error = sqrt(error_gpu);
  //CkPrintf("[Iter %d, Chare %d-%d-%d] error: %f\n", iter, xc, yc, zc, error);
  if(error < derefine_cutoff && thisIndex.getDepth() > min_depth) return COARSEN;
  else if(error > refine_cutoff && thisIndex.getDepth() < max_depth) return REFINE;
  else return STAY;
#else // USE_HAPI
  // create callback
  CkArrayIndexOctIndex myIndex = CkArrayIndexOctIndex(thisIndex);
  CkCallback *cb = new CkCallback(CkIndex_MeshBlock::gotErrorFromGPU(), myIndex, thisProxy);

  // offload
  invokeDecisionKernel(decisionStream, u, h_error, d_error, d_u, mesh_manager_local->d_delu, mesh_manager_local->d_delua, refine_filter, dx, dy, dz, block_x, cb);

  return STAY; // return dummy value
#endif // USE_HAPI
#endif // USE_GPU
}

void MeshBlock::resetMeshRestructureData(){
  decision = INV;
  parentHasAlreadyMadeDecision=false;
  hasReceivedParentDecision=false;

  uncleDecisions.clear();

  for (std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(),
       iend = neighbors.end(); it != iend; ++it) {
    it->second.resetDecision();
  }

  for (int i = 0; i < NUM_CHILDREN; ++i)
    child_decision[i] = INV;
}

void MeshBlock::makeGranularityDecisionAndCommunicate(){
  if(isLeaf) {//run this on leaf nodes
//#ifndef USE_HAPI
#if 1 // TODO Don't use HAPI version because it sometimes results in different refinement decisions
    Decision newDecision = (decision != REFINE) ? std::max(decision, getGranularityDecision()) : decision;
    VB(logfile << thisIndex.getIndexString().c_str() << " decision = " << newDecision << std::endl;);
    updateDecisionState(1, newDecision);
#else
    // With HAPI the decision is made asynchronously, so we can't decide yet
    getGranularityDecision();
#endif
  }
  else if(isGrandParent() && !parentHasAlreadyMadeDecision) informParent(gen_iter,-1, INV, 1);
}

/***** PHASE1 FUNCTIONS****/
void MeshBlock::updateDecisionState(int cascade_length, Decision newDecision) {
  cascade_length++;

  assert(isLeaf);

  if (decision == newDecision) return;

  assert((decision != REFINE || !isRefined) && "Re-refining?!");

  decision = newDecision;
  if (decision == COARSEN) return; // Don't communicate the 'default' decision

  for (int i=0; i<NUM_NEIGHBORS; i++) {
    OctIndex QI = thisIndex.getNeighbor(i);

    if (neighbors.count(QI)) {
      Neighbor & N = neighbors[QI];
      if (!N.isRefined()) {
        // isFriend
        VB(logfile << "sending exchangePhase1Msg, decision = " << decision << ", receiver = " << QI.getIndexString() << std::endl;);
        thisProxy(QI).exchangePhase1Msg(gen_iter, getSourceDirection(i), -1, decision, cascade_length);
      } else {
        // isNephew
        //Get Corresponding Children of the neighbor
        std::vector<OctIndex> children;
        getChildrenInDir(QI, getSourceDirection(i), children);
        // XXX: -2 means "we're your uncle"
        for (std::vector<OctIndex>::iterator I = children.begin(),
             E = children.end(); I != E; ++I){
        VB(logfile << "sending exchangePhase1Msg, decision = " << decision << ", receiver = " << I->getIndexString() << std::endl;);
        thisProxy(*I).exchangePhase1Msg(gen_iter, getSourceDirection(i), -2, decision, cascade_length);
        }
      }
    } else {
      // Does not exist, talk to uncle
      VB(logfile << "sending exchangePhase1Msg, decision = " << decision << ", receiver = " << QI.getParent().getIndexString() << std::endl;);
      thisProxy(QI.getParent()).exchangePhase1Msg(gen_iter, getSourceDirection(i), thisIndex.getOctant(), decision, cascade_length);
    }
  }

  if(parent != thisIndex){
    VB(logfile << "sending informParent, decision = " << decision << std::endl;);
    thisProxy(parent).informParent(gen_iter, thisIndex.getOctant(), decision, cascade_length);
  }
}

// Will be called from two contexts:
//  a) If the parent is a grandparent and also have children that are leaves
//  b) when a child sends REFINE/STAY message to the parent
void MeshBlock::processChildDecision(int childNum, Decision dec, int cascade_length) {
  VB(logfile << "recvd informParent, childNum " << childNum << std::endl;);
  if(childNum >= 0) child_decision[childNum]=dec;

  if(dec==REFINE) child_isRefined[childNum]=true;
  
  if(parentHasAlreadyMadeDecision == false){
    parentHasAlreadyMadeDecision = true;
    FOR_EACH_CHILD
      if(i!=childNum && !child_isRefined[i]){
        VB(logfile << "sending message to child " << thisIndex.getChild(i).getIndexString().c_str() << ", miterattions " << gen_iter << std::endl;);
        thisProxy(thisIndex.getChild(i)).recvParentDecision(gen_iter, cascade_length);
      }
    END_FOR
    if(parent!=thisIndex)
      thisProxy(parent).informParent(gen_iter, thisIndex.getOctant(), STAY, cascade_length);
  }
}

void MeshBlock::processParentDecision(int cascade_length) {
  hasReceivedParentDecision = true;
  Decision newDecision = std::max(STAY, decision);
  if(isLeaf) updateDecisionState(cascade_length, newDecision);
}

void MeshBlock::processPhase1Msg(int dir, int quadrant, Decision remoteDecision, int cascade_length) {
  VB(logfile << "isLeaf: " << isLeaf << std::endl;);
  Decision newDecision = decision;
  OctIndex QI = thisIndex.getNeighbor(dir);
  VB(logfile << "received exchangePhase1Msg, dir " << dir << ", quadrant " << quadrant << ", idx " << QI.getIndexString() << std::endl;);
  VB(logfile << QI.getIndexString() << " decision: " << remoteDecision << std::endl;);
  if(!isLeaf)
    logfile << "only leaves should exchange remesh messages" << std::endl;
  assert(isLeaf && "Only leaves should exchange remesh messages");

  if (quadrant == -2)
    assert(!neighbors.count(QI) && "Uncle in our neighbor map");
  else{
    VB(if(neighbors.count(QI)==0){
      logfile << thisIndex.getIndexString().c_str() << " didn't knew that " << QI.getIndexString().c_str() << " existed" << std::endl;

      for(std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        logfile << it->first.getIndexString().c_str() << std::endl;
    });
    assert(neighbors.count(QI) && "Received message from neighbor we didn't know existed!");
  }

  // Ignore decisions sent from our uncle
  if (!neighbors.count(QI)) {
    VB(logfile << "received message from uncle" << std::endl;);
    OctIndex uncleIndex = QI.getParent();
    if (uncleDecisions.count(uncleIndex)) {
      uncleDecisions[uncleIndex] = std::max(uncleDecisions[uncleIndex], remoteDecision);
    } else {
      uncleDecisions[uncleIndex] = remoteDecision;
    }
    VB(logfile << "done processing message: decision = " << decision << std::endl;);
    return;
  }

  Neighbor &N = neighbors[QI];

  if(quadrant == -1)
    assert(!N.isRefined());
  else
    assert(N.isRefined());
  remoteDecision = N.setDecision(remoteDecision, quadrant);

  if(!N.isRefined()) {
    if(remoteDecision == REFINE) {
      newDecision = std::max(STAY, decision);
    }
  } else {
    newDecision = std::max(decision, remoteDecision);
  }

  updateDecisionState(cascade_length, newDecision);
}

/**** PHASE2 FUNCTIONS ****/
void MeshBlock::doPhase2(){
  VB(logfile << "in doPhase2, iteration = " << iter << " decision = " << decision << std::endl;);
  if(decision == COARSEN){//send data to the parent
    std::vector<float> child_u;
    if(!inInitialMeshGenerationPhase){
      child_u.resize((block_y*block_x*block_z)/8);
      for(int i=1; i<= block_x; i+=2)
        for(int j=1; j<=block_y; j+=2)
          for(int k=1; k<=block_z; k+=2)
            child_u[index_c(i/2, j/2, k/2)] = downSample(u, i, j, k);
    }
    //CkPrintf("[Iter %d, Depth %d, Chare %d-%d-%d] coarsening\n", iter, thisIndex.getDepth(), xc, yc, zc);
    VB(logfile << "coarsening .. sending data to parent " << gen_iter << std::endl;);
    thisProxy(parent).recvChildData(gen_iter, thisIndex.getOctant(), myt, mydt, gen_iter, iter, child_u, neighbors, uncleDecisions);
    thisProxy[thisIndex].ckDestroy();
    //thisProxy.doneInserting();
    return;
  }
  else if(decision == REFINE) refine();

  updateMeshState();
  resetMeshRestructureData();
  //iterate();
}

void MeshBlock::updateMeshState(){
  //Update the Status of Your Neighbors, need to be done only if you are going to stay in that position
  if(isLeaf && decision == STAY){

    std::vector<OctIndex> ToRemove;
    for (std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(),
         iend = neighbors.end(); it != iend; ++it) {
      Neighbor &N = it->second;
      if (!N.isRefined()) {
        switch(N.getDecision()) {
          case REFINE:      N.setRefined(true); break;
          case STAY:        N.setRefined(false); break;
          case COARSEN:     ToRemove.push_back(it->first); break;
          //default:          CkAbort("nbr_decision not set");
        }
      }
      else {
        // TODO: Only check nephews sharing a surface with us
        // For now, abuse the fact that default is COARSEN
        // and only nephews adjacent to us send messages,
        // so the rest of the data is invalid anyway.
        Decision D = COARSEN;
        for (int i = 0; i < NUM_CHILDREN; ++i)
          D = std::max(D, N.getDecision(i));
        switch (D) {
          case STAY:      N.setRefined(true); break;
          case COARSEN:   N.setRefined(false); break;
          case REFINE:    CkAbort("unacceptable decision");
          //default:        CkAbort("nbr_decision not set");
        }
      }
    }

    for (int i = 0; i < NUM_NEIGHBORS; ++i) {
      OctIndex QI = thisIndex.getNeighbor(i);
      OctIndex uncleIndex = QI.getParent();
      if (uncleDecisions.count(uncleIndex)) {
        assert(!neighbors.count(uncleIndex));
        assert(!neighbors.count(QI));

        Decision d = uncleDecisions[uncleIndex];
        if (d == REFINE) {
          neighbors[QI] = Neighbor(i);
        }
      }
    }

    // Remove the elements marked, avoid iterator badness
    for (std::vector<OctIndex>::iterator I = ToRemove.begin(),
         E = ToRemove.end(); I != E; ++I) {
      VB(logfile << "erasing " << I->getIndexString().c_str() << std::endl;);
      assert(neighbors.count(*I));
      neighbors.erase(*I);
    }
  }
  else if(decision == REFINE){
    isRefined=true;
  }
  else if(isRefined && !isGrandParent() && !parentHasAlreadyMadeDecision){
    VB(logfile << thisIndex.getIndexString().c_str() << " setting isRefined to false "
          << isRefined << " " << isGrandParent() << " " << parentHasAlreadyMadeDecision << std::endl;);
    isRefined = false;
    //set up the sdag for receving child data
    thisProxy[thisIndex].wait4ChildData();
    return;
  }

  if(isGrandParent()) {
    FOR_EACH_CHILD
      if(child_decision[i]==INV && child_isRefined[i])//did not receive any message
        child_isRefined[i]=false;
    END_FOR
  }
  iterate();
}

void MeshBlock::recvChildData(int childNum, float myt, float mydt,
                              int gen_iter, int iter, std::vector<float> child_u,
                              std::map<OctIndex, Neighbor> childNeighbors,
                              std::map<OctIndex, Decision> childUncleDecisions){
  VB(logfile << "recvd data from child: " << childNum << std::endl;);
  this->myt = myt;
  this->mydt = mydt;
  this->iter = iter;
  this->gen_iter = gen_iter;

  int st_i, end_i, st_j, end_j, st_k, end_k;
  getOctantRange(childNum, st_i, end_i, st_j, end_j, st_k, end_k);

  int ctr=0; float rsq;

  // If a leaf is coarsening into us, we should know
  // that we are no longer refined.  Check this now.
  if(isRefined) VB(logfile << "leaves coarsened into refined node" << std::endl;);
  assert(!isRefined && "Leaves coarsened into refined node?");

  if(inInitialMeshGenerationPhase) applyInitialCondition();
  else
    for (int k=st_k; k<=end_k; k++)
      for(int j=st_j; j<=end_j; j++){
        for(int i=st_i; i<=end_i; i++)
          u[index(i,j,k)]=child_u[ctr++];
      }

  //Update the Status of Your childNeighbors based on Data Sent from the Children

  int c1, c2, c3;
  switch(childNum){
    case 0: c1=LEFT;  c2=DOWN; c3=BACKWARD; break;
    case 1: c1=LEFT;  c2=DOWN; c3=FORWARD;  break;
    case 2: c1=LEFT;  c2=UP;   c3=BACKWARD; break;
    case 3: c1=LEFT;  c2=UP;   c3=FORWARD;  break;
    case 4: c1=RIGHT; c2=DOWN; c3=BACKWARD; break;
    case 5: c1=RIGHT; c2=DOWN; c3=FORWARD;  break;
    case 6: c1=RIGHT; c2=UP;   c3=BACKWARD; break;
    case 7: c1=RIGHT; c2=UP;   c3=FORWARD;  break;
  }

  this->iter  = iter;
  this->myt         = myt;
  this->mydt        = mydt;
  setNbrStateUponCoarsening(c1, childNum, childNeighbors, childUncleDecisions);
  setNbrStateUponCoarsening(c2, childNum, childNeighbors, childUncleDecisions);
  setNbrStateUponCoarsening(c3, childNum, childNeighbors, childUncleDecisions);

  resetMeshRestructureData();
}

inline void MeshBlock::setNbrStateUponCoarsening(int dir, int childNum, std::map<OctIndex, Neighbor> & childNeighbors, std::map<OctIndex, Decision> & childUncleDecisions) {
  OctIndex QIParent = thisIndex.getNeighbor(dir);
  OctIndex QIChild = thisIndex.getChild(childNum).getNeighbor(dir);

  if (childNeighbors.count(QIChild)) {
    Neighbor parentNeighbor = Neighbor(dir);
    Neighbor childNeighbor = childNeighbors[QIChild];

    if(childNeighbor.isRefined()) {
      parentNeighbor.setRefined(true);
    } else {
      switch(childNeighbor.getDecision()) {
        case COARSEN: parentNeighbor.setRefined(false); break;
        case STAY:    parentNeighbor.setRefined(true); break;
      }
    }

    // assert(!neighbors.count(QIParent));
    neighbors[QIParent] = parentNeighbor;
  } else if(neighbors.count(QIParent)) {
    Neighbor& N = neighbors[QIParent];

    // Okay, last we heard we do have a neighbor
    // in this direction, and because the first branch above
    // wasn't taken we know that this neighbor doesn't have a child.
    // Question is: does it even exist?
    // Look to our child (that is just not merging into us) for
    // if it received a decision (non-COARSEN), and act accordingly.
    // If it didn't receive a decision, the neighbor's decision
    // must have been to coarsen, and remove it.
    Decision D = COARSEN;
    if (childUncleDecisions.count(QIParent)) {
      D = childUncleDecisions[QIParent];
    }
    switch(D) {
      case REFINE: N.setRefined(true); break;
      case STAY: N.setRefined(false); break;
      case COARSEN: neighbors.erase(QIParent);
    }
  } else if (childUncleDecisions.count(QIParent)) {
    // Our child has information about our neighbor,
    // but we DON'T.  This is okay, this can happen
    // if *our* uncle refined in the same iteration
    // we did, as we only update our 'neighbor state'
    // if we're STAY.  Anyway, look at what our
    // neighbor decided to do and update accordingly:
    switch(childUncleDecisions[QIParent]) {
      case STAY:
        neighbors[QIParent] = Neighbor(dir);
        break;
      case REFINE:
        neighbors[QIParent] = Neighbor(dir);
        neighbors[QIParent].setRefined(true);
        break;
      case COARSEN: break; // We already don't know about it
    }
  }
}

void MeshBlock::interpolate(float *u, std::vector<float>& refined_u, int xstart, int xend, int ystart, int yend, int zstart, int zend){
  float sx_l, sx_r, sy_u, sy_d, sz_b, sz_f;
  for(int i=xstart, m=1; i<=xend; i++, m++){
    for(int j=ystart, n=1; j<=yend; j++, n++){
      for(int k=zstart, o=1; k<=zend; k++, o++){
        sx_l = (-u[index(i-1,j,k)]+u[index(i,j,k)]) / 4;
        sx_r = (-u[index(i,j,k)]+u[index(i+1,j,k)]) / 4;

        sy_d = (-u[index(i,j-1,k)]+u[index(i,j,k)]) / 4;
        sy_u = (-u[index(i,j,k)]+u[index(i,j+1,k)]) / 4;

        sz_b = (-u[index(i,j,k-1)]+u[index(i,j,k)]) / 4;
        sz_f = (-u[index(i,j,k)]+u[index(i,j,k+1)]) / 4;

        refined_u[index_l(2*(m-1),   2*(n-1)  , 2*(o-1))] = u[index(i,j,k)] - sx_l - sy_d - sz_b;
        refined_u[index_l(2*(m-1)+1, 2*(n-1)  , 2*(o-1))] = u[index(i,j,k)] + sx_r - sy_d - sz_b;
        refined_u[index_l(2*(m-1),   2*(n-1)+1, 2*(o-1))] = u[index(i,j,k)] - sx_l + sy_u - sz_b;
        refined_u[index_l(2*(m-1)+1, 2*(n-1)+1, 2*(o-1))] = u[index(i,j,k)] + sx_r + sy_u - sz_b;
        refined_u[index_l(2*(m-1),   2*(n-1)  , 2*(o-1)+1)] = u[index(i,j,k)] - sx_l - sy_d + sz_f;
        refined_u[index_l(2*(m-1)+1, 2*(n-1)  , 2*(o-1)+1)] = u[index(i,j,k)] + sx_r - sy_d + sz_f;
        refined_u[index_l(2*(m-1),   2*(n-1)+1, 2*(o-1)+1)] = u[index(i,j,k)] - sx_l + sy_u + sz_f;
        refined_u[index_l(2*(m-1)+1, 2*(n-1)+1, 2*(o-1)+1)] = u[index(i,j,k)] + sx_r + sy_u + sz_f;
      }
    }
  }
}

void MeshBlock::refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, int zstart, int zend, float x_min, float y_min, float z_min) {
  OctIndex child = thisIndex.getChild(sChild);

  std::vector<float> refined_u;
  if(!inInitialMeshGenerationPhase){
    refined_u.resize(block_x*block_y*block_z);
    interpolate(u, refined_u, xstart, xend, ystart, yend, zstart, zend);
  }
  VB(logfile << thisIndex.getIndexString().c_str() << " isRefined = " << isRefined << std::endl; 
  logfile << thisIndex.getIndexString().c_str() << " inserting " << child.getIndexString().c_str() << std::endl;);

  // Creation of new chares due to refinement
  //CkPrintf("[Iter %d, Depth %d, Chare %d-%d-%d] refining\n", iter, thisIndex.getDepth(), xc, yc, zc);
  thisProxy(child).insert(dx/2, dy/2, dz/2, myt, mydt, x_min, y_min, z_min, gen_iter, iter, refined_u, neighbors);
}

void MeshBlock::refine(){

  for (unsigned c = 0; c < NUM_CHILDREN; ++c) {
    int cx_min, cx_max, cy_min, cy_max, cz_min, cz_max;
    getOctantRange(c, cx_min, cx_max, cy_min, cy_max, cz_min, cz_max);

    float cxx = x_min, cyy = y_min, czz = z_min;
    if (cx_max == block_x)
      cxx += (nx*dx)/2;
    if (cy_max == block_y)
      cyy += (ny*dy)/2;
    if (cz_max == block_z)
      czz += (nz*dz)/2;

    refineChild(c, cx_min, cx_max, cy_min, cy_max, cz_min, cz_max, cxx, cyy, czz);
    //thisProxy.doneInserting();
  }
  allFree();
}

bool MeshBlock::isGrandParent() {
  bool ret = false;
  for (int i = 0; i < NUM_CHILDREN; ++i) ret = ret || child_isRefined[i];
  return ret;
}


void MeshBlock::startLdb(){
  UserSetLBLoad();
  /*if(isRoot()){ 
    ckout << "at sync" << endl;
  }*/
  //thisArray->remoteDoneInserting();
  //ckout << thisIndex.getIndexString().c_str() << " in startLdb " << isLeaf << " " << CkMyPe() << endl;
  AtSync();
  //ckout << "at sync" << endl;
  //if (isRoot())
  VB(logfile << "starting ldb" << std::endl;);
  VB(logfile << "memory usage: " << CmiMemoryUsage() << std::endl;);
  if(isRoot()){
    //LBDatabaseObj()->StartLB();
    ckout << "ldb start time: " << CmiWallTimer() << endl;
  }
}

void MeshBlock::ResumeFromSync() {
  if(isRoot()) ckout << "ldb end time: " << CmiWallTimer() << endl;
  VB(logfile << "resuming from load balancing" << std::endl;);
  //ckout <<  thisIndex.getIndexString().c_str() << " " << isLeaf << endl;
  if(isLeaf)
    mesh_manager_local->incrementWorkUnitCount(iter);
  startPhase2(gen_iter);
}

bool MeshBlock::isRoot() {
  return thisIndex.nbits == min_depth * numDims && thisIndex.bitVector == 0;
}

void MeshBlock::printData() {
  int cntr=0;
  for(int xIndex = 1; xIndex <= block_x; xIndex++){
    for(int yIndex = 1; yIndex <= block_y; yIndex++) {
      for(int zIndex = 1; zIndex <= block_z; zIndex++) {
        logfile << x[xIndex] << " " << y[yIndex] << " " << z[zIndex] << " " << u[index(xIndex,yIndex,zIndex)] << " " << cntr++ << std::endl;
      }
    }
  }

  //thisProxy[parent].donePrinting();
}

#include "Advection.def.h"
