#ifndef HEADERS_H_
#define HEADERS_H_

#include <limits>
#include <fstream>
#include <cmath>
#include <cstring>
#include <queue>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

#include "charm++.h"
#include "pup_stl.h"

using namespace std;

#ifdef LOGGER
#define VB(x)  x
#else
#define VB(x)
#endif

#define NUM_NEIGHBORS (2 * numDims)
#define NUM_CHILDREN (1 << numDims)

#define FOR_EACH_ZONE for (int i = 1; i <= block_x; i++) \
                        for (int j = 1; j <= block_y; j++) \
                          for (int k = 1; k <= block_z; k++) {
#define FOR_EACH_NEIGHBOR for (int i = 0; i < NUM_NEIGHBORS; i++) {
#define FOR_EACH_CHILD for (int i = 0; i < NUM_CHILDREN; i++) {
#define END_FOR }

#define isLeaf !isRefined

extern int block_x, block_y, block_z;

static inline int index2d(int a, int b, int width) { return b * width + a; }
static inline int index_xz(int x, int z) { return z * block_x + x; }
static inline int index_xy(int x, int y) { return y * block_x + x; }
static inline int index_yz(int y, int z) { return z * block_y + y; }

static inline int index(int i, int j, int k) {
  return (k * (block_y + 2) + j) * (block_x + 2) + i;
}

static inline int index_l(int i, int j, int k) {
  return (k * block_y + j) * block_x + i;
}

static inline int index_c(int i, int j, int k) {
  return (k * (block_y / 2) + j) * (block_x / 2) + i;
}

static inline bool liesin(int x, int st, int end) {
  return x <= end && x >= st;
}

const int numBitsForStorage = 8 * sizeof(unsigned int);
const int numDims = 3;
const int numDims2 = numDims * numDims;

enum Dir { UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3, FORWARD = 4, BACKWARD = 5 };
PUPbytes(Dir);

enum Decision { INV = -1, COARSEN = 0, STAY = 1, REFINE = 2 };
PUPbytes(Decision);

#include "OctIndex.h"
#include "Advection.h"
#include "Main.h"

#endif // HEADERS_H_
