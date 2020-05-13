#include "Headers.h"
#include <assert.h>

extern CProxy_Main mainProxy;
extern CProxy_Advection qtree;

extern int array_height, array_width, array_depth;
extern int block_width, block_height, block_depth;
extern int num_chare_rows, num_chare_cols, num_chare_Zs;
extern int min_depth, max_depth;

extern int max_iters, refine_freq, lb_freq;

extern bool verbose;

extern float vx, vy, vz;
extern float apx, anx, apy, any, apz, anz;
extern float x_ctr, y_ctr, z_ctr, radius;
extern float tmax, t, dt, cfl;

#define inInitialMeshGenerationPhase (meshGenIterations <= max_depth)

float refine_filter = 0.01;
float refine_cutoff= 0.2, derefine_cutoff = 0.05;

CProxy_MeshManager mesh_manager;

#ifdef USE_GPU
extern void memHostAlloc(void** ptr, size_t size);
extern void memHostFree(void* ptr);
extern void memDeviceAlloc(void** ptr, size_t size);
extern void memDeviceFree(void* ptr);
extern void createStream(cudaStream_t* stream_ptr);
extern void destroyStream(cudaStream_t stream);

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
inline static void getChildrenInDir(OctIndex myIndex, int dir, std::vector<OctIndex>& children)
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

  for (std::vector<int>::iterator iter = ids.begin(), end = ids.end(); iter != end; ++iter)
    children.push_back(myIndex.getChild(*iter));
  assert(children.size() == 4);
}

inline static void setFirstHalf(int& min, int& max) { max /= 2; }
inline static void setSecondHalf(int& min, int& max) { min = (max / 2) + 1; }

inline static void populateQuadrant(bool bit1, bool bit2, int& min1, int& max1, int& min2, int& max2)
{
  if (bit1)  setSecondHalf(min1, max1);
  else       setFirstHalf(min1, max1);

  if (bit2)  setSecondHalf(min2, max2);
  else       setFirstHalf(min2, max2);
}

enum {
  X_MASK = 1 << 2,
  Y_MASK = 1 << 1,
  Z_MASK = 1 << 0
};

inline static void populateRanges(int dir, int octant, int &x_min, int &x_max,
                           int& y_min, int& y_max, int& z_min, int& z_max)
{
  x_min = 0; x_max = block_width-1;
  y_min = 0; y_max = block_height-1;
  z_min = 0; z_max = block_depth-1;

  switch (dir) {
    case UP:
      y_min = y_max = block_height+1;
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
      x_min = x_max = block_width+1;
      if (octant >= 0) {
        populateQuadrant(octant & Y_MASK, octant & Z_MASK, y_min, y_max, z_min, z_max);
      }
      break;
    case FORWARD:
      z_min = z_max = block_depth+1;
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
  x_min = 1; x_max = block_width;
  y_min = 1; y_max = block_height;
  z_min = 1; z_max = block_depth;

  if (octant & X_MASK)  x_min = block_width/2+1;
  else                  x_max = block_width/2;
  if (octant & Y_MASK)  y_min = block_height/2+1;
  else                  y_max = block_height/2;
  if (octant & Z_MASK)  z_min = block_depth/2+1;
  else                  z_max = block_depth/2;
}

MeshManager::MeshManager() : workUnitCount(0), compute_time_sum(0.0),
  compute_time_cnt(0), decision_time_sum(0.0), decision_time_cnt(0)
{
  // delu and delua are 4D arrays
#ifdef USE_GPU
  size_t delu_size = sizeof(float)*numDims*(block_width+2)*(block_height+2)*(block_depth+2);
  memDeviceAlloc((void**)&d_delu, delu_size);
  memDeviceAlloc((void**)&d_delua, delu_size);
#endif
  delu = new float***[numDims];
  delua = new float***[numDims];

  for (int d = 0; d < numDims; ++d) {
    delu[d] = new float**[block_width+2];
    delua[d] = new float**[block_width+2];
    for (int i = 0; i < block_width+2; ++i) {
      delu[d][i] = new float*[block_height+2];
      delua[d][i] = new float*[block_height+2];
      for (int j = 0; j < block_height+2; ++j) {
        delu[d][i][j] = new float[block_depth+2];
        delua[d][i][j] = new float[block_depth+2];
      }
    }
  }
}

MeshManager::MeshManager(CkMigrateMessage* m) : CBase_MeshManager(m)
{
#ifdef USE_GPU
  size_t delu_size = sizeof(float)*numDims*(block_width+2)*(block_height+2)*(block_depth+2);
  memDeviceAlloc((void**)&d_delu, delu_size);
  memDeviceAlloc((void**)&d_delua, delu_size);
#endif
  delu = new float***[numDims];
  delua = new float***[numDims];

  for (int d = 0; d < numDims; ++d) {
    delu[d] = new float**[block_width+2];
    delua[d] = new float**[block_width+2];
    for (int i = 0; i < block_width+2; ++i) {
      delu[d][i] = new float*[block_height+2];
      delua[d][i] = new float*[block_height+2];
      for (int j = 0; j < block_height+2; ++j) {
        delu[d][i][j] = new float[block_depth+2];
        delua[d][i][j] = new float[block_depth+2];
      }
    }
  }
}

MeshManager::~MeshManager()
{
#ifdef USE_GPU
  memDeviceFree(d_delu);
  memDeviceFree(d_delua);
#endif
}

void MeshManager::incrementWorkUnitCount(int iterations)
{
  workUnitCount++;
  workUnits[iterations]++;

  if (iterations % lb_freq == 0 || iterations % lb_freq == lb_freq-1) {
    // This is either a load balancing iteration or the one right before it
    minLoad[iterations] += 1;
    maxLoad[iterations] += 1;
    avgLoad[iterations] += 1;
  }
}

void MeshManager::reduceWorkUnits()
{
  CkCallback cb(CkReductionTarget(Main,totalWorkUnits), mainProxy);
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

void Advection::prepareData4Exchange(){
  imsg=0;
  ghostReceived.clear();

  for (std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(),
       iend = neighbors.end(); it != iend; ++it) {
      it->second.setDataSent(false);
  }

  // Cache boundary data in contiguous blocks

  // YZ surfaces
  for(int j=1; j <= block_height; ++j)
    for(int k=1; k <= block_depth; ++k) {
      left_surface[index_yz(j-1,k-1)] = u[index(1,j,k)];
      right_surface[index_yz(j-1,k-1)] = u[index(block_width,j,k)];
    }

  // XZ Surfaces
  for(int i=1; i <= block_width; ++i)
    for(int k=1; k <= block_depth; ++k) {
      top_surface[index_xz(i-1,k-1)] = u[index(i, block_height, k)];
      bottom_surface[index_xz(i-1,k-1)] = u[index(i, 1, k)];
    }

  // XY surfaces
  for(int i=1; i <= block_width; ++i)
    for(int j=1; j <= block_height; ++j) {
      forward_surface[index_xy(i-1,j-1)] = u[index(i, j, block_depth)];
      backward_surface[index_xy(i-1,j-1)] = u[index(i, j, 1)];
    }
}

void Advection::applyInitialCondition(){
  VB(logfile << "applying initial condition" << std::endl;);
  float rcub;
  for(int i=0; i<block_width+2; i++)
    for(int j=0; j<block_height+2; j++)
      for(int k=0; k<block_depth+2; k++){
        rcub = (x[i] - x_ctr)*(x[i]-x_ctr) +
               (y[j] - y_ctr)*(y[j]-y_ctr) +
               (z[k] - z_ctr)*(z[k]-z_ctr);
        u[index(i, j, k)] = (rcub<=radius*radius) ? 2:1;
        //VB(logfile << x[i] << " " << y[j] << " " << z[k] << ": " << rcub << " " << radius*radius << " " << u[index(i, j, k)] << std::endl;);
      }
}

void Advection::mem_allocate(float* &p, int size){
  p = new float[size];
}

void Advection::mem_allocate_all(){
#ifdef USE_GPU
  memHostAlloc((void**)&u, (block_width+2)*(block_height+2)*(block_depth+2)*sizeof(float));
  memHostAlloc((void**)&u2, (block_width+2)*(block_height+2)*(block_depth+2)*sizeof(float));
  memHostAlloc((void**)&u3, (block_width+2)*(block_height+2)*(block_depth+2)*sizeof(float));
  memHostAlloc((void**)&h_error, sizeof(float));
  memDeviceAlloc((void**)&d_u, sizeof(float)*(block_width+2)*(block_height+2)*(block_depth+2));
  memDeviceAlloc((void**)&d_u2, sizeof(float)*(block_width+2)*(block_height+2)*(block_depth+2));
  memDeviceAlloc((void**)&d_u3, sizeof(float)*(block_width+2)*(block_height+2)*(block_depth+2));
  memDeviceAlloc((void**)&d_error, sizeof(float));

  createStream(&computeStream);
  createStream(&decisionStream);
#else
  mem_allocate(u, (block_width+2)*(block_height+2)*(block_depth+2));
  mem_allocate(u2, (block_width+2)*(block_height+2)*(block_depth+2));
  mem_allocate(u3, (block_width+2)*(block_height+2)*(block_depth+2));
#endif

  mem_allocate(x, block_width+2);
  mem_allocate(y, block_height+2);
  mem_allocate(z, block_depth+2);

  mem_allocate(left_surface, block_height*block_depth);
  mem_allocate(right_surface, block_height*block_depth);
  mem_allocate(top_surface, block_width*block_depth);
  mem_allocate(bottom_surface, block_width*block_depth);
  mem_allocate(forward_surface, block_width*block_height);
  mem_allocate(backward_surface, block_width*block_height);
}

void Advection::mem_deallocate_all(){
#ifdef USE_GPU
  memHostFree(u);
  memHostFree(u2);
  memHostFree(u3);
  memHostFree(h_error);
  memDeviceFree(d_u);
  memDeviceFree(d_u2);
  memDeviceFree(d_u3);
  memDeviceFree(d_error);

  destroyStream(computeStream);
  destroyStream(decisionStream);
#else
  delete [] u;
  delete [] u2;
  delete [] u3;
#endif

  delete [] x;
  delete [] y;
  delete [] z;

  delete [] left_surface;
  delete [] right_surface;
  delete [] top_surface;
  delete [] bottom_surface;
  delete [] forward_surface;
  delete [] backward_surface;
}

// Constructor used to create initial chares
Advection::Advection(float x_min, float x_max, float y_min, float y_max,
                     float z_min, float z_max)
{
  mesh_manager_local = mesh_manager.ckLocalBranch();

  thisIndex.getCoordinates(xc, yc, zc);
  dx = (x_max - x_min)/float(array_width);
  dy = (y_max - y_min)/float(array_height);
  dz = (z_max - z_min)/float(array_depth);

  nx = array_width/(num_chare_cols);
  ny = array_height/(num_chare_rows);
  nz = array_depth/(num_chare_Zs);

  myt = t;
  mydt = dt;

  this->x_min = xc*nx*dx;
  this->y_min = yc*ny*dy;
  this->z_min = zc*nz*dz;

  iterations=0;
  meshGenIterations=0;

  // Allocate all necessary buffers
  mem_allocate_all();

  initializeRestofTheData();
}

void Advection::initializeRestofTheData(){
  usesAutoMeasure = false;
  usesAtSync = true;
  remeshStartTime = 0;
  VB(logfile.open(string("log/"+thisIndex.getIndexString()+"log").c_str()););
  if(inInitialMeshGenerationPhase){
    for(int i=0; i<block_width+2; i++)
      x[i] = x_min + float(i)*dx - 0.5*dx;


    for(int i=0; i<block_height+2; i++)
      y[i] = y_min + float(i)*dy - 0.5*dy;

    for(int i=0; i<block_depth+2; i++)
      z[i] = z_min + float(i)*dz - 0.5*dz;

    applyInitialCondition();
  }

  FOR_EACH_CHILD
    child_isRefined[i] = false;
  END_FOR
  this->isRefined = false;

  FOR_EACH_NEIGHBOR
    VB(logfile << "neighbor in dir " << i << " is " << thisIndex.getNeighbor(i).getIndexString().c_str() << std::endl;);
    neighbors[thisIndex.getNeighbor(i)] = Neighbor(i);
  END_FOR

  parent = (thisIndex.getDepth()>min_depth) ? thisIndex.getParent() : thisIndex;
  resetMeshRestructureData();
  nChildDataRecvd = 0;
  phase1Over = false;
}

//added for array migration - see how 2D arrays can be packed
void Advection::pup(PUP::er &p){
  p|isRefined;
  p|depth;

  PUParray(p, child_isRefined, NUM_CHILDREN);

  p|ghostReceived;

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
  p|imsg;
  p|nChildDataRecvd;
  p|phase1Over;

  if(p.isUnpacking()){
    mem_allocate_all(); // Needed as there is no migration constructor
    //resetMeshRestructureData();
  }

  PUParray(p, u, (block_width+2)*(block_height+2)*(block_depth+2));
  PUParray(p, x, block_width+2);
  PUParray(p, y, block_height+2);
  PUParray(p, z, block_depth+2);

  p|iterations;
  p|meshGenIterations;
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
  p|nChildDataRecvd;
  p|phase1Over;
  p|amr3d_i;
  p|ichild;
}

Advection::~Advection(){
  if (isLeaf) {
    mem_deallocate_all(); // FIXME: Shouldn't this be called for ALL chares?
  }
}

inline float downSample(float* u, int x, int y, int z) {
  return (
    u[index(x, y, z)]   + u[index(x+1, y, z)] +
    u[index(x, y+1, z)] + u[index(x+1, y+1, z)] +
    u[index(x, y, z+1)]   + u[index(x+1, y, z+1)] +
    u[index(x, y+1, z+1)] + u[index(x+1, y+1, z+1)]) / 8.0;
}

class surface_iterator {
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
  surface_iterator(float *u,
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

  inline surface_iterator& operator++() {
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

float* Advection::getGhostBuffer(int dir) {
  switch (dir) {
    case UP:    return top_surface;
    case DOWN:  return bottom_surface;
    case LEFT:  return left_surface;
    case RIGHT: return right_surface;
    case FORWARD: return forward_surface;
    case BACKWARD: return backward_surface;
  }
}

int Advection::getGhostCount(int dir) {
  switch (dir) {
  case UP:    case DOWN: return block_width*block_depth;
  case RIGHT: case LEFT: return block_height*block_depth;
  case FORWARD: case BACKWARD: return block_height*block_width;
  //default: CkAbort("Asking for an unknown boundary's size");
  }
}

void Advection::sendGhost(int dir){
  int count = getGhostCount(dir);
  VB(logfile << "ghost count in direction " << dir << " = " << count << std::endl;);
  float* boundary;

  OctIndex QI = thisIndex.getNeighbor(dir);
  std::map<OctIndex, Neighbor>::iterator I = neighbors.find(QI);

  if (I == neighbors.end()) {
    VB(logfile << "neighbor is an uncle" << std::endl;);
    // Uncle case, neighbor doesn't exist in this direction (at this level)
    OctIndex receiver = QI.getParent();

    // TODO: Remove use of subdirections
    boundary = getGhostBuffer(dir);
    count /= 4;

    int x_min, x_max, y_min, y_max, z_min, z_max;
    populateRanges(dir, -1, x_min, x_max, y_min, y_max, z_min, z_max);

    int dx = (x_min == x_max) ? 0 : 2;
    int dy = (y_min == y_max) ? 0 : 2;
    int dz = (z_min == z_max) ? 0 : 2;
    if(x_min!=x_max)x_min++;
    else if(x_min==0) x_min = x_max = 1;
    else if(x_min==block_width+1) x_min = x_max = block_width-1;

    if(y_min!=y_max)y_min++;
    else if(y_min==0) y_min = y_max = 1;
    else if(y_min==block_width+1) y_min = y_max = block_height-1;

    if(z_min!=z_max)z_min++;
    else if(z_min==0) z_min = z_max = 1;
    else if(z_min==block_width+1) z_min = z_max = block_depth-1;

    surface_iterator iter(u, x_min, x_max, dx,
                             y_min, y_max, dy,
                             z_min, z_max, dz);
    int k;
    for (k = 0; !iter.isDone(); ++iter, ++k)
      boundary[k] = downSample(u, iter.getX(), iter.getY(), iter.getZ());
    assert(k == count);

    thisProxy(receiver).receiveGhosts(iterations, getSourceDirection(dir), thisIndex.getOctant(), count, boundary);
  } else if (!I->second.isRefined()) {
    // Friend
    VB(logfile << "neighbor " << QI.getIndexString() << " is a friend" << std::endl;);
    boundary = getGhostBuffer(dir);
    thisProxy(QI).receiveGhosts(iterations, getSourceDirection(dir), -1, count, boundary);
  }
}

void Advection::process(int iteration, int dir, int quadrant, int size, float gh[]){
  VB(logfile << "received ghost from direction " << dir << ", octant " << quadrant << std::endl;);
  bool fromNephew = (quadrant >= 0);

  OctIndex QI = thisIndex.getNeighbor(dir);
  if (fromNephew) QI = QI.getChild(quadrant);
  //assert(!ghostReceived.count(QI));
  ghostReceived.insert(QI);

  imsg += (fromNephew) ? 0.25:1;

  int x_min, x_max, y_min, y_max, z_min, z_max;
  populateRanges(dir, quadrant, x_min, x_max, y_min, y_max, z_min, z_max);
  VB(logfile << "after populate ranges " << " (iteration : " << iteration << ") " << dir << " " << quadrant << " " \
             << x_min << " " << x_max << " " << y_min << " " << y_max << " " \
             << z_min << " " << z_max << std::endl;);
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

  surface_iterator iter(u, x_min, x_max, dx,
                           y_min, y_max, dy,
                           z_min, z_max, dz);

  for(int i=0; i<size; ++i, ++iter){
      VB(
        if(y[iter.getY()]==0.46875) 
            //logfile << std::setw(5) << gh[i] << " ";
            logfile << "setting " << x[iter.getX()] << ", " << y[iter.getY()] << ", " << z[iter.getZ()] << " = " << gh[i] << std::endl;
       ); 
      *iter = gh[i];
  }
  VB(logfile << std::endl;);
  assert(iter.isDone());
}

void Advection::sendReadyData(){
  //check if data can be sent to any of the refined neighbors
  //If the neighbors are at the same level or do not exist at all 
  //data will be sent in begin_iteration function and need 
  //not be sent here
  for (std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(),
       iend = neighbors.end(); it != iend; ++it) {
    Neighbor &N = it->second;
    if(N.isRefined() && !N.isDataSent()) {
      std::vector<OctIndex> children;
      getChildrenInDir(thisIndex, N.getDir(), children);

      std::vector<OctIndex> allNeighbors;
      for (std::vector<OctIndex>::iterator I = children.begin(),
           E = children.end(); I != E; ++I) {
        std::vector<OctIndex> neighbors = I->getNeighbors();
        allNeighbors.insert(allNeighbors.end(),
                            neighbors.begin(), neighbors.end());
      }

      // TODO: change to iterator
      bool receivedAll = true;
      for(int i = 0; i < allNeighbors.size(); i++) {
        if(allNeighbors[i].getParent() != thisIndex) { // external neighbor
          if(!ghostReceived.count(allNeighbors[i]) && !ghostReceived.count(allNeighbors[i].getParent())) {
            receivedAll = false;
            break;
          }
        }
      }

      if(receivedAll) {
        interpolateAndSend(N.getDir());
        N.setDataSent(true);
      }
    }
  }
}

int Advection::getSourceDirection(int NBR) {
  switch (NBR) {
  case UP:    return DOWN;
  case DOWN:  return UP;
  case LEFT:  return RIGHT;
  case RIGHT: return LEFT;
  case FORWARD: return BACKWARD;
  case BACKWARD: return FORWARD;
  }
}
    
void Advection::interpolateAndSend(int dir) {
  int uncledir = getSourceDirection(dir);
  std::vector<OctIndex> children;
  getChildrenInDir(thisIndex.getNeighbor(dir), uncledir, children);
  for (std::vector<OctIndex>::iterator I = children.begin(),
       E = children.end(); I != E; ++I)
    interpolateAndSendToNephew(uncledir, *I);
}

void Advection::interpolateAndSendToNephew(int uncledir, OctIndex QI) {
  float *boundary;
  int count = getGhostCount(uncledir);

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
    columncount = block_height;
    dx = 0;
    y_min++; y_max++;
    z_min++; z_max++;
    if(x_min==0){
        x_min = x_max = 1;
    }else{
        x_min = x_max = block_width;
    }

  }
  if (y_min == y_max) {
    a = sx;
    b = sz;
    c = (y_max == 0) ? &sy_d : &sy_u;
    columncount = block_width;
    dy = 0;
    x_min++; x_max++;
    z_min++; z_max++;
    if(y_min==0){
        y_min = y_max = 1;
    }else{
        y_min = y_max = block_height;
    }
  }
  if (z_min == z_max) {
    a = sx;
    b = sy;
    c = (z_max == 0) ? &sz_b : &sz_f;
    columncount = block_width;
    dz = 0;
    x_min++; x_max++;
    y_min++; y_max++;
    if(z_min==0){
        z_min = z_max = 1;
    }else{
        z_min = z_max = block_depth;
    }
  }

  surface_iterator in(u, x_min, x_max, dx,
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

  thisProxy(QI).receiveGhosts(iterations, uncledir, -2, count, boundary);
}

void Advection::computeDone() {
#if defined(USE_GPU) && defined(USE_HAPI)
  double time_dur = CkWallTimer() - compute_start_time;
  mesh_manager_local->addComputeTime(time_dur);

  iterate();
#endif
}

void Advection::compute(){
  //if(iterations==1){
#ifdef LOGGER
    char logfilename[100];
    sprintf(logfilename, "out/snap_%s_%d.vtk", thisIndex.getIndexString().c_str(), iterations);
    logfile.open(logfilename);
    logfile.precision(8);
    //logfile << "Variables=\"X\",\"Y\",\"Z\",\"radius\",\"nodeid\"\n";
    //logfile << "Zone I = 16, J = 16, K = 16, F = POINT\n";
    printData();
    logfile.close();
    /*int dims[] = {block_height+1, block_width+1, block_depth+1};
    int vardims[] = {1};
    int centering[] = {0};
    const char *varnames[]  = {"viscosity"};
    //varnames = new char*[1];
    //varnames[0] = "viscosity";
    
    float *vars[1];
    vars[0] = new float[block_width*block_height*block_depth];
    for(int i=0; i<block_width; i++)
        for(int j=0; j<block_height; j++)
            for(int k=0; k<block_depth; k++){
                vars[0][(k*block_height+j)*block_width+i] = u[index(i+1,j+1,k+1)];
                //ckout << x[i+1] << " " << y[j+1] << " " << z[k+1] << " " << u[index(i+1,j+1,k+1)] << endl;
            }
    float *xn = new float[block_width+1];
    float *yn = new float[block_height+1];
    float *zn = new float[block_depth+1];
    for(int i=0; i<block_width+1; i++) {xn[i]=x_min+dx*i;}
    for(int i=0; i<block_height+1; i++){yn[i]=y_min+dy*i;}
    for(int i=0; i<block_depth+1; i++) {zn[i]=z_min+dz*i;}
    //float *vars[] = {u};
    write_rectilinear_mesh(logfilename, 0, dims, xn, yn, zn, 1, vardims, centering, varnames, vars);
    delete [] vars[0]; delete [] xn; delete [] yn; delete [] zn;
    for(int k=0; k<=block_depth+1; k++){
        for(int j=0; j<=block_height+1; j++){
            for(int i=0; i<=block_width+1; i++){
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
  memcpy(u2, u, sizeof(float)*(block_width+2)*(block_height+2)*(block_depth+2));
  memcpy(u3, u, sizeof(float)*(block_width+2)*(block_height+2)*(block_depth+2));
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
  invokeComputeKernel(computeStream, u, d_u, d_u2, d_u3, dx, dy, dz, dt, apx, apy, apz, anx, any, anz, block_width, NULL);
#else
  // create callback
  CkArrayIndexOctIndex myIndex = CkArrayIndexOctIndex(thisIndex);
  CkCallback *cb = new CkCallback(CkIndex_Advection::computeDone(), myIndex, thisProxy);

  invokeComputeKernel(computeStream, u, d_u, d_u2, d_u3, dx, dy, dz, dt, apx, apy, apz, anx, any, anz, block_width, cb);
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

void Advection::gotErrorFromGPU() {
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

Decision Advection::getGranularityDecision(){
  float delx = 0.5/dx;
  float dely = 0.5/dy;
  float delz = 0.5/dz;
  float error=0;

  decision_start_time = CkWallTimer();

//#ifndef USE_GPU
#if 1 // FIXME Do not use the GPU kernels for now, the error values are slightly different
  /********** CPU CODE **********/
  for(int i=1; i <= block_width; i++){
    for(int j=1; j<=block_height; j++){
      for(int k=1; k<=block_depth; k++){
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

  int istart=2, iend=block_width-1,
      jstart=2, jend=block_height-1,
      kstart=2, kend=block_depth-1;
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
  //CkPrintf("[Iter %d, Chare %d-%d-%d] error: %f\n", iterations, xc, yc, zc, error);
  if(error < derefine_cutoff && thisIndex.getDepth() > min_depth) return COARSEN;
  else if(error > refine_cutoff && thisIndex.getDepth() < max_depth) return REFINE;
  else return STAY;
#else
  /********** GPU CODE **********/
//#ifndef USE_HAPI
#if 1 // FIXME Don't use HAPI version because it sometimes results in different refinement decisions
  // execute GPU kernel
  float error_gpu = invokeDecisionKernel(decisionStream, u, h_error, d_error, d_u, mesh_manager_local->d_delu, mesh_manager_local->d_delua, refine_filter, dx, dy, dz, block_width, NULL);

  double decision_time = CkWallTimer() - decision_start_time;
  mesh_manager_local->addDecisionTime(decision_time);

  error = sqrt(error_gpu);
  //CkPrintf("[Iter %d, Chare %d-%d-%d] error: %f\n", iterations, xc, yc, zc, error);
  if(error < derefine_cutoff && thisIndex.getDepth() > min_depth) return COARSEN;
  else if(error > refine_cutoff && thisIndex.getDepth() < max_depth) return REFINE;
  else return STAY;
#else // USE_HAPI
  // create callback
  CkArrayIndexOctIndex myIndex = CkArrayIndexOctIndex(thisIndex);
  CkCallback *cb = new CkCallback(CkIndex_Advection::gotErrorFromGPU(), myIndex, thisProxy);

  // offload
  invokeDecisionKernel(decisionStream, u, h_error, d_error, d_u, mesh_manager_local->d_delu, mesh_manager_local->d_delua, refine_filter, dx, dy, dz, block_width, cb);

  return STAY; // return dummy value
#endif // USE_HAPI
#endif // USE_GPU
}

void Advection::resetMeshRestructureData(){
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

void Advection::makeGranularityDecisionAndCommunicate(){
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
  else if(isGrandParent() && !parentHasAlreadyMadeDecision) informParent(meshGenIterations,-1, INV, 1);
}

/***** PHASE1 FUNCTIONS****/
void Advection::updateDecisionState(int cascade_length, Decision newDecision) {
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
        thisProxy(QI).exchangePhase1Msg(meshGenIterations, getSourceDirection(i), -1, decision, cascade_length);
      } else {
        // isNephew
        //Get Corresponding Children of the neighbor
        std::vector<OctIndex> children;
        getChildrenInDir(QI, getSourceDirection(i), children);
        // XXX: -2 means "we're your uncle"
        for (std::vector<OctIndex>::iterator I = children.begin(),
             E = children.end(); I != E; ++I){
        VB(logfile << "sending exchangePhase1Msg, decision = " << decision << ", receiver = " << I->getIndexString() << std::endl;);
        thisProxy(*I).exchangePhase1Msg(meshGenIterations, getSourceDirection(i), -2, decision, cascade_length);
        }
      }
    } else {
      // Does not exist, talk to uncle
      VB(logfile << "sending exchangePhase1Msg, decision = " << decision << ", receiver = " << QI.getParent().getIndexString() << std::endl;);
      thisProxy(QI.getParent()).exchangePhase1Msg(meshGenIterations, getSourceDirection(i), thisIndex.getOctant(), decision, cascade_length);
    }
  }

  if(parent != thisIndex){
    VB(logfile << "sending informParent, decision = " << decision << std::endl;);
    thisProxy(parent).informParent(meshGenIterations, thisIndex.getOctant(), decision, cascade_length);
  }
}

// Will be called from two contexts:
//  a) If the parent is a grandparent and also have children that are leaves
//  b) when a child sends REFINE/STAY message to the parent
void Advection::processChildDecision(int childNum, Decision dec, int cascade_length) {
  VB(logfile << "recvd informParent, childNum " << childNum << std::endl;);
  if(childNum >= 0) child_decision[childNum]=dec;

  if(dec==REFINE) child_isRefined[childNum]=true;
  
  if(parentHasAlreadyMadeDecision == false){
    parentHasAlreadyMadeDecision = true;
    FOR_EACH_CHILD
      if(i!=childNum && !child_isRefined[i]){
        VB(logfile << "sending message to child " << thisIndex.getChild(i).getIndexString().c_str() << ", miterattions " << meshGenIterations << std::endl;);
        thisProxy(thisIndex.getChild(i)).recvParentDecision(meshGenIterations, cascade_length);
      }
    END_FOR
    if(parent!=thisIndex)
      thisProxy(parent).informParent(meshGenIterations, thisIndex.getOctant(), STAY, cascade_length);
  }
}

void Advection::processParentDecision(int cascade_length) {
  hasReceivedParentDecision = true;
  Decision newDecision = std::max(STAY, decision);
  if(isLeaf) updateDecisionState(cascade_length, newDecision);
}

void Advection::processPhase1Msg(int dir, int quadrant, Decision remoteDecision, int cascade_length) {
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
void Advection::doPhase2(){
  VB(logfile << "in doPhase2, iteration = " << iterations << " decision = " << decision << std::endl;);
  if(decision == COARSEN){//send data to the parent
    std::vector<float> child_u;
    if(!inInitialMeshGenerationPhase){
      child_u.resize((block_height*block_width*block_depth)/8);
      for(int i=1; i<= block_width; i+=2)
        for(int j=1; j<=block_height; j+=2)
          for(int k=1; k<=block_depth; k+=2)
            child_u[index_c(i/2, j/2, k/2)] = downSample(u, i, j, k);
    }
    //CkPrintf("[Iter %d, Depth %d, Chare %d-%d-%d] coarsening\n", iterations, thisIndex.getDepth(), xc, yc, zc);
    VB(logfile << "coarsening .. sending data to parent " << meshGenIterations << std::endl;);
    thisProxy(parent).recvChildData(meshGenIterations, thisIndex.getOctant(), myt, mydt, meshGenIterations, iterations, child_u, neighbors, uncleDecisions);
    thisProxy[thisIndex].ckDestroy();
    //thisProxy.doneInserting();
    return;
  }
  else if(decision == REFINE) refine();

  updateMeshState();
  resetMeshRestructureData();
  //iterate();
}

void Advection::updateMeshState(){
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

void Advection::recvChildData(int childNum, float myt, float mydt,
                              int meshGenIterations, int iterations, std::vector<float> child_u,
                              std::map<OctIndex, Neighbor> childNeighbors,
                              std::map<OctIndex, Decision> childUncleDecisions){
  VB(logfile << "recvd data from child: " << childNum << std::endl;);
  this->myt = myt;
  this->mydt = mydt;
  this->iterations = iterations;
  this->meshGenIterations = meshGenIterations;

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

  this->iterations  = iterations;
  this->myt         = myt;
  this->mydt        = mydt;
  setNbrStateUponCoarsening(c1, childNum, childNeighbors, childUncleDecisions);
  setNbrStateUponCoarsening(c2, childNum, childNeighbors, childUncleDecisions);
  setNbrStateUponCoarsening(c3, childNum, childNeighbors, childUncleDecisions);

  resetMeshRestructureData();
  //if (nChildDataRecvd++ == NUM_CHILDREN){
  //  nChildDataRecvd = 0;
  //  iterate();
  //}
}

inline void Advection::setNbrStateUponCoarsening(int dir, int childNum, std::map<OctIndex, Neighbor> & childNeighbors, std::map<OctIndex, Decision> & childUncleDecisions) {
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

void Advection::interpolate(float *u, std::vector<float>& refined_u, int xstart, int xend, int ystart, int yend, int zstart, int zend){
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

void Advection::refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, int zstart, int zend, float x_min, float y_min, float z_min) {
  OctIndex child = thisIndex.getChild(sChild);

  std::vector<float> refined_u;
  if(!inInitialMeshGenerationPhase){
    refined_u.resize(block_width*block_height*block_depth);
    interpolate(u, refined_u, xstart, xend, ystart, yend, zstart, zend);
  }
  VB(logfile << thisIndex.getIndexString().c_str() << " isRefined = " << isRefined << std::endl; 
  logfile << thisIndex.getIndexString().c_str() << " inserting " << child.getIndexString().c_str() << std::endl;);

  // Creation of new chares due to refinement
  //CkPrintf("[Iter %d, Depth %d, Chare %d-%d-%d] refining\n", iterations, thisIndex.getDepth(), xc, yc, zc);
  thisProxy(child).insert(dx/2, dy/2, dz/2, myt, mydt, x_min, y_min, z_min, meshGenIterations, iterations, refined_u, neighbors);
}

void Advection::refine(){

  for (unsigned c = 0; c < NUM_CHILDREN; ++c) {
    int cx_min, cx_max, cy_min, cy_max, cz_min, cz_max;
    getOctantRange(c, cx_min, cx_max, cy_min, cy_max, cz_min, cz_max);

    float cxx = x_min, cyy = y_min, czz = z_min;
    if (cx_max == block_width)
      cxx += (nx*dx)/2;
    if (cy_max == block_height)
      cyy += (ny*dy)/2;
    if (cz_max == block_depth)
      czz += (nz*dz)/2;

    refineChild(c, cx_min, cx_max, cy_min, cy_max, cz_min, cz_max, cxx, cyy, czz);
    //thisProxy.doneInserting();
  }
  mem_deallocate_all();
}

bool Advection::isGrandParent() {
  bool ret = false;
  for (int i = 0; i < NUM_CHILDREN; ++i) ret = ret || child_isRefined[i];
  return ret;
}

// Constructor used to create children chares on refinement
Advection::Advection(float dx, float dy, float dz,
                     float myt, float mydt, float x_min, float y_min, float z_min,
                     int meshGenIterations, int iterations, std::vector<float> refined_u, std::map<OctIndex, Neighbor> parentNeighbors)
{
  mesh_manager_local = mesh_manager.ckLocalBranch();

  this->dx = dx;
  this->dy = dy;
  this->dz = dz;

  this->myt = myt;
  this->mydt = mydt;

  this->x_min = x_min;
  this->y_min = y_min;
  this->z_min = z_min;

  nx = array_width/(num_chare_cols);
  ny = array_height/(num_chare_rows);
  nz = array_depth/(num_chare_Zs);

  thisIndex.getCoordinates(xc, yc, zc);
  this->meshGenIterations = meshGenIterations;
  this->iterations = iterations;

  // Allocate all necessary buffers
  mem_allocate_all();

  initializeRestofTheData();

  nChildDataRecvd=0;
  std::map<OctIndex, Neighbor>::iterator it, iend;
  for (it = neighbors.begin(), iend = neighbors.end(); it != iend; ) {
    const OctIndex &neighborOctIndex = it->first;
    Neighbor &neighbor = it->second;
    bool shouldDelete = false;
    if (neighborOctIndex.getParent() != thisIndex.getParent()) {
      std::map<OctIndex, Neighbor>::iterator parentIt = parentNeighbors.find(neighborOctIndex.getParent());
      if (parentIt == parentNeighbors.end()) {
        shouldDelete = true;
      } else {
        Neighbor parentNeighbor = parentIt->second;
        if (parentNeighbor.isRefined()) {
          switch (parentNeighbor.getDecision(neighborOctIndex.getOctant())) {
            case COARSEN: shouldDelete = true; break;
            //case STAY: break;
            case REFINE: neighbor.setRefined(true); break;
          }
        } else {
          switch (parentNeighbor.getDecision()) {
            //case REFINE: break;
            case STAY: shouldDelete = true; break;
            case COARSEN: VB(logfile << "this neighbor cannot derefine" << std::endl;); 
                        CkAbort("this neighbor cannot derefine");
          }
        }
      }
    }

    // Advance iterator, removing this quadindex/neighbor pair if requested.
    if (shouldDelete)
      neighbors.erase(it++);
    else
      ++it;
  }

  assert(neighbors.size() >= numDims);
  unsigned sameParent = 0;
  for (std::map<OctIndex, Neighbor>::iterator it = neighbors.begin(),
       iend = neighbors.end(); it != iend; ++it) {
    if (it->first.getParent() == thisIndex.getParent())
      ++sameParent;
  }
  assert(sameParent == 3);

  if(!inInitialMeshGenerationPhase){
    int ctr=0;
    for(int k=1; k<=block_depth; k++)
      for(int j=1; j<=block_height; j++)
        for(int i=1; i<=block_width; i++)
          u[index(i,j,k)]=refined_u[ctr++];
  }
  iterate();
}

void Advection::startLdb(){
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

void Advection::ResumeFromSync() {
  if(isRoot()) ckout << "ldb end time: " << CmiWallTimer() << endl;
  VB(logfile << "resuming from load balancing" << std::endl;);
  //ckout <<  thisIndex.getIndexString().c_str() << " " << isLeaf << endl;
  if(isLeaf)
    mesh_manager_local->incrementWorkUnitCount(iterations);
  startPhase2(meshGenIterations);
}

bool Advection::isRoot(){
  return thisIndex.nbits == min_depth * numDims && thisIndex.bitVector == 0;
}

void Advection::printData() {
  int cntr=0;
  for(int xIndex = 1; xIndex <= block_width; xIndex++){
    for(int yIndex = 1; yIndex <= block_height; yIndex++) {
      for(int zIndex = 1; zIndex <= block_depth; zIndex++) {
        logfile << x[xIndex] << " " << y[yIndex] << " " << z[zIndex] << " " << u[index(xIndex,yIndex,zIndex)] << " " << cntr++ << std::endl;
      }
    }
  }

  //thisProxy[parent].donePrinting();
}

#include "Advection.def.h"
