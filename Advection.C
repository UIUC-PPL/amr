#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include<limits>
using namespace std;

#include "charm++.h"
#include "liveViz.h"
#include "pup_stl.h"
#include "trace-projections.h"

#include "Constants.h"
#include "QuadIndex.h"
#include "Main.decl.h"
#include "Main.h"

extern CProxy_Main mainProxy;

extern int array_height;
extern int array_width;

extern int num_chare_rows;
extern int num_chare_cols;

extern int block_width, block_height;

extern int min_depth, max_depth;

#define wrap_x(a) (((a)+num_chare_rows)%num_chare_rows)
#define wrap_y(a) (((a)+num_chare_cols)%num_chare_cols)


extern int nframe;
extern double xctr, yctr, radius;
extern double v;
extern double ap, an;
extern double tmax, t, dt, cfl;
extern int max_iterations, refine_frequency;

#define index(i,j)  (int)((j)*(block_width+2) + i)

const int ndim=2;
const int ndim2 = 4; //ndim*ndim

double ***delu, ***delua;
double delu2[ndim2], delu3[ndim2], delu4[ndim2];
double delx, dely, dely_f;
double refine_filter = 0.01;
double refine_cutoff=0.8, derefine_cutoff=0.2;

#include "Advection.decl.h"
#include "Advection.h"
PerProcessorChare::PerProcessorChare(){
  delu = new double**[ndim];
  delua = new double**[ndim];

  delu[0] = new double*[block_height+2];
  delu[1] = new double*[block_height+2];
  delua[0] = new double*[block_height+2];
  delua[1] = new double*[block_height+2];


  for(int i=0; i<block_height+2; i++){
      delu[0][i] = new double[block_width+2];
      delu[1][i] = new double[block_width+2];
      delua[0][i] = new double[block_width+2];
      delua[1][i] = new double[block_width+2];
  }
}

InitRefineMsg::InitRefineMsg(bool isInMeshGenerationPhase, double dx, double dy, 
                             double myt, double mydt, double xmin, double  ymin, 
                             int meshGenIterations_, int iterations_, vector<double>& refined_u, 
                             bool *nbr_exists,
                             bool *nbr_isRefined, DECISION *nbr_decision) {
  this->isInMeshGenerationPhase = isInMeshGenerationPhase;
  this->meshGenIterations = meshGenIterations_;
  this->dx = dx;
  this->dy = dy;
  this->myt = myt;
  this->mydt = mydt;
  this->xmin = xmin;
  this->ymin = ymin;
  this->iterations = iterations_;
  if (!isInMeshGenerationPhase) 
    memcpy(this->refined_u, &refined_u[0], sizeof(double)*block_height*block_width);
  memcpy(this->parent_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
  memcpy(this->parent_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
  memcpy(this->parent_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);
}

void Advection::applyInitialCondition(){
  double rsq;
  for(int i=0; i<block_width+2; i++){
    for(int j=0; j<block_height+2; j++){
      rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[j]-yctr);
        if(rsq <= radius*radius)
          u[index(i, block_height+1-j)] = 2;
        else u[index(i, block_height+1-j)] = 1;
    }
  }
  return;

  double *xarray = x;
  double *yarray = y;

  double x, y, rx, ry;
  double t = 2;
  for(int i = 0; i < block_width + 2; i++){
    for(int j = 0; j < block_height + 2; j++){
      x = xarray[i]*100; y = yarray[j]*100;
      u[index(i, block_height+1-j)] = 1;
      if (!(y >= 42 && y < 58))
          continue;
      if (x >= 8 && x < 16){//its the C
        //get corresponding y1max, y1min and y2ax, y2min
        double y1max, y1min, y2min, y2max;
        y1max = 50 + sqrt(64 - pow((x - 16), 2.));
        y1min = y1max - t;
        y2min = 50 - sqrt(64 - pow((x - 16), 2.));
        y2max = y2min + t;
        if ((y < y1max && y >= y1min) || (y >= y2min && y < y2max))
          u[index(i, block_height+1-j)] = 2;
      }else if(x >= 20 && x < 28){//its H
        rx = x - 20; ry = y - 42;
        if(rx <= 2 || rx >= 6 || (ry >= 6 && ry < 10))
          u[index(i, block_height+1-j)] = 2;
      }else if(x >= 32 && x < 40){//its A
         rx = x - 32; ry = y - 42;
         if ((ry >= 8 && ry < 8 + t) || (ry >= 14))
           u[index(i, block_height+1-j)] = 2;
         if ((rx < t) || (rx >= 8 - t))
           u[index(i, block_height+1-j)] = 2;
      }else if(x >= 44 && x < 52){//its R
         rx = x - 44; ry = y - 42;
         if (rx < t || (rx >= 6 && ry >= 8))
           u[index(i, block_height+1-j)] = 2;
         if ((ry >= 8 && ry <= 10) || (ry >= 14))
           u[index(i, block_height+1-j)] = 2;
         if (rx >= 2){
            double ymax = -1*(4.0/3.0)*rx + (32./3.);
            double ymin = ymax - t;
            if (ry >= ymin && ry < ymax)
              u[index(i, block_height+1-j)] = 2;
         }
      }else if(x >= 56 && x < 64){//its M
        rx = x - 56; ry = y - 42;
        if ((rx < t) || (rx >= 8 - t))
           u[index(i, block_height+1-j)] = 2;
        if (ry >= 12)
           u[index(i, block_height+1-j)] = 2;
        if (rx >= 3 && rx < 5 && ry >= 8)
           u[index(i, block_height+1-j)] = 2;
      }else if(x >= 68 && x < 76){//its +
        rx = x - 69; ry = y - 42;
        if ((ry >= 5 && ry <= 10) || (rx >= 2 && rx <= 4))
          u[index(i, block_height+1-j)] = 2;
      }else if(x >= 80 && x < 88){//its +
        rx = x - 81; ry = y - 42;
        if ((ry >= 5 && ry <= 10) || (rx >= 2 && rx <= 4))
          u[index(i, block_height+1-j)] = 2;
      }
    }
  }
}

void Advection::mem_allocate(double* &p, int size){
  p = new double[size];
}


void Advection::mem_allocate_all(){
  VB(logFile << "Allocating Memory for " << thisIndex.getIndexString() << std::endl;);
  mem_allocate(u, (block_width+2)*(block_height+2));
  mem_allocate(u2, (block_width+2)*(block_height+2));
  mem_allocate(u3, (block_width+2)*(block_height+2));

  mem_allocate(x, block_width+2);
  mem_allocate(y, block_height+2);

  mem_allocate(left_edge, block_height);
  mem_allocate(right_edge, block_height);

  mem_allocate(top_edge, block_width);
  mem_allocate(bottom_edge, block_width);
}

Advection::Advection(double xmin, double xmax, double ymin, double ymax)
  /*: AdvTerm(thisProxy, thisIndex, true)*/
{
  //ckout << thisIndex.getIndexString().c_str() << " created" << endl;
//Constructor for the Initial Grid Zones
  __sdag_init();

  usesAtSync = CmiTrue;

  has_terminated=false;
  hasReset = false;
  shouldDestroy = false;
  char fname[100];
  sprintf(fname, "log/%s.log", thisIndex.getIndexString().c_str());
  VB(logFile.open(fname););
    
    //srand(thisIndex.getQuadI() + atoi(thisIndex.getIndexString()));

  this->isRefined = false;
    
  for(int dir=UP; dir<=RIGHT; ++dir)
    nbr[dir] = thisIndex.getNeighbor(dir);
  /*    
        if(thisIndex.nbits!=0)
        parent = thisIndex.getParent();
        else parent = thisIndex;*/

  /*Nodes constructed here will not have any parent as this 
    is the coarsest level of the Mesh*/
  parent = thisIndex;
  isGrandParent = false;
  /*Lets set the nieghbors status if this is the root node */
  for(int i=0; i<NUM_NEIGHBORS; i++){
    nbr_exists[i]=true;
    nbr_isRefined[i]=false;
    nbr_dataSent[i]=false;
  }
    
  int xc, yc;
  thisIndex.getCoordinates(xc, yc);

  dx = (xmax - xmin)/double(array_height);
  dy = (ymax - ymin)/double(array_width);

  nx = array_height/(num_chare_cols);
  ny = array_width/(num_chare_rows);

  myt = t;
  mydt = dt;

  this->xmin = xc*nx*dx;
  this->ymin = yc*ny*dy;

  for (int i = 0; i < NUM_CHILDREN; ++i)
    child_isRefined[i] = false;

  VB(logFile << "xmin: " << this->xmin << ", ymin: " << this->ymin << std::endl;);

    advection();
}

void Advection::printState(){
  QuadIndex qindex = thisIndex;
  VB(logFile << "Printing Status of Node " << qindex.getIndexString() << std::endl;);
    VB(logFile << "isRefined: " << isRefined << std::endl;
       logFile << "nbr_exists: ";);
#ifdef LOGGER
    for(int j=0; j<NUM_NEIGHBORS; j++)
      logFile << nbr_exists[j] << ", ";
  logFile << std::endl << "nbr_isRefined: ";
  for(int j=0; j<NUM_NEIGHBORS; j++)
    logFile << nbr_isRefined[j] << ", ";
  logFile << std::endl;
#endif
}

void Advection::advection(){

  //logFile << "In Advection: " << std::endl;
  myt = t;
  mem_allocate_all();
  iterations=0;
  meshGenIterations=0;

  thisIndex.getCoordinates(xc, yc);

  for(int i=0; i<block_width+2; i++){
    x[i] = xmin + double(i)*dx - 0.5*dx;
    //logFile << x[i] << std::endl;
  }
  //logFile << std::endl;

  for(int i=0; i<block_height+2; i++){
    y[i] = ymin + double(i)*dy - 0.5*dy;
    // logFile << y[i] << std::endl;
  }
  //logFile << std::endl;

  double rsq;
    
  //logFile << "In Adfvection2" << std::endl;
  VB(logFile << "xctr: " << xctr << ", yctr: " << yctr << ", radius: " << radius << std::endl;);
  applyInitialCondition();
#if 1
#ifdef LOGGER
  for(int i=0; i<block_height; i++){
    for(int j=0; j<block_width; j++)
      logFile << u[index(j+1,i+1)] << "\t";
    logFile << std::endl;
  }
  logFile << std::endl;
#endif 
  //CkExit();
#endif
  char fname[100];
  sprintf(fname, "out/out_%s_%d", thisIndex.getIndexString().c_str(), iterations);
  VB(outFile.open(fname););

#ifdef LOGGER

    //outFile << "coordinates: " << xc << ", " << yc << std::endl;
    //outFile << dx << ", " << dy << std::endl;
    for(int i=1; i<=block_width; i++){
      for(int j=1; j<=block_height; j++){
        //outFile << xmin << ", " << double(xc*block_width + i) << std::endl;
        outFile << xmin + (double(i))*dx - 0.5*dx << " "\
                << ymin + (double(j))*dy - 0.5*dy << " "\
                << u[index(i,block_height+1-j)] << std::endl;
      }
    }
  outFile.flush();
  outFile.close();
#endif
}

//added for array migration - see how 2D arrays can be packed
void Advection::pup(PUP::er &p){
  VB(logFile << "In PUP" << std::endl;);
    CBase_Advection::pup(p);
  __sdag_pup(p);

  p|isRefined;
  p|isGrandParent;
  p|depth;

  for(int i=0; i<NUM_CHILDREN; i++)
    p|child_isRefined[i];

  for(int i=0; i<NUM_NEIGHBORS; i++){
    p|nbr_exists[i];
    p|nbr_isRefined[i];
  }

  for(int i=0; i<3*NUM_NEIGHBORS; i++){
    p|nbr_dataSent[i];
  }

  p|hasReceived;
  p|hasReset;
  p|decision;
  p|parentHasAlreadyMadeDecision;
  p|hasReceivedParentDecision;
  p|hasCommunicatedSTAY;
  p|hasCommunicatedREFINE;
  p|hasAllocatedMemory;

  p|parent;
  for(int i=0; i<NUM_NEIGHBORS; i++)
      p|nbr[i];

  p|xc;
  p|yc;
  p|imsg;

  if(p.isUnpacking())
    mem_allocate_all();
    
  for(int i=0; i<(block_width+2)*(block_height+2); i++){
    p|u[i];
    //p|u2[i];
    //p|u3[i];
  }

  for (int i=0; i<block_width+2; i++){
    p|x[i];
    //p|top_edge[i];
    //p|bottom_edge[i];
  }

  for (int i=0; i<block_height+2; i++){
    p|y[i];
    //p|left_edge[i];
    //p|right_edge[i];
  }
 
  p|iterations;
  p|meshGenIterations;
  p|shouldDestroy;
  p|up;
  p|un;
  p|myt;
  p|mydt;

  p|dx; p|dy; p|nx; p|ny; p|xmin; p|xmax; p|ymin; p|ymax;
    
}
    
Advection::~Advection(){
  //CkPrintf("%s destructor %d\n", thisIndex.getIndexString().c_str(), iterations);
  VB(logFile << "In Destructor" << std::endl;);
  delete [] u;
  delete [] u2;
  delete [] u3;

  delete [] x;
  delete [] y;
  delete [] left_edge;
  delete [] right_edge;

  delete [] top_edge;
  delete [] bottom_edge;
}

inline double downSample(double* u, int x, int y) {
  return (u[index(x, y)]   + u[index(x+1, y)] +
          u[index(x, y+1)] + u[index(x+1, y+1)]) / 4;
}

struct boundary_iterator {
  double *u;
  int x, y, dx, dy;

  inline double left()  {    return u[index(x-1, y)];  }
  inline double right() {    return u[index(x+1, y)];  }
  inline double up()    {    return u[index(x, y+1)];  }
  inline double down()  {    return u[index(x, y-1)];  }
  inline double& operator*() {
    return u[index(x, y)];
  }
  inline boundary_iterator& operator++() {
    x += dx;
    y += dy;
    return *this;
  }
  inline bool operator!=(const boundary_iterator& rhs) {
    CkAssert(u == rhs.u);
    return x != rhs.x || y != rhs.y;
  }
  inline bool operator==(const boundary_iterator& rhs) {
    CkAssert(u == rhs.u);
    return x == rhs.x && y == rhs.y;
  }
  inline boundary_iterator operator+(int offset) {
    boundary_iterator ret = *this;
    ret.x += dx * offset;
    ret.y += dy * offset;
    return ret;
  }
  boundary_iterator(double* u_, int x_, int dx_, int y_, int dy_)
    : u(u_), x(x_), y(y_), dx(dx_), dy(dy_)
  { }
  boundary_iterator() {}
};

double* Advection::getGhostBuffer(int dir) {
  switch (dir) {
  case UP:    case UP_LEFT:   case UP_RIGHT:   return top_edge;
  case DOWN:  case DOWN_LEFT: case DOWN_RIGHT: return bottom_edge;
  case LEFT:  case LEFT_UP:   case LEFT_DOWN:  return left_edge;
  case RIGHT: case RIGHT_UP:  case RIGHT_DOWN: return right_edge;
  default: CkAbort("Asking for an unknown boundary");
  }
}

int Advection::getGhostCount(int dir) {
  switch (dir) {
  case UP:    case UP_LEFT:   case UP_RIGHT:
  case DOWN:  case DOWN_LEFT: case DOWN_RIGHT: return block_width;
  case LEFT:  case LEFT_UP:   case LEFT_DOWN:
  case RIGHT: case RIGHT_UP:  case RIGHT_DOWN: return block_height;
  default: CkAbort("Asking for an unknown boundary's size");
  }
}

bool Advection::sendGhost(int dir, bool which=0){//which = 0 - message sent during remeshing, which=1 - during normal course of iteration
  int count = getGhostCount(dir);
  double* boundary;

  if(nbr_exists[dir] && !nbr_isRefined[dir]){
    int val = rand();
    VB(logFile << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Neighbor " << nbr[dir].getIndexString() << ", iteration " << iterations << ", " << SENDER_DIR[dir] << ", " << val << std::endl;);

    switch (dir) {
    case LEFT:  boundary = left_edge;                  break;
    case RIGHT: boundary = right_edge;                 break;
    case UP:    boundary = &u[index(1, 1)];            break;
    case DOWN:  boundary = &u[index(1, block_height)]; break;
    default:
      CkPrintf("noop send\n");
      return false;
    };
    thisProxy(nbr[dir]).receiveGhosts(iterations, SENDER_DIR[dir], count, boundary,
                                        thisIndex, val);
    lastSent = nbr[dir];
    return true;
  }

  if(!nbr_exists[dir]) {
    QuadIndex receiver = thisIndex.getNeighbor(dir).getParent();
    int sender_direction = map_nbr(thisIndex.getQuadI(), dir);
    boundary = getGhostBuffer(dir);
    VB(logFile << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Uncle: " << receiver.getIndexString() << ", iteration " << iterations << endl;);
    count /= 2;

    boundary_iterator begin, end;
    switch(dir) {
    case UP:    begin = boundary_iterator(u, 1,               2, 1,                0); break;
    case DOWN:  begin = boundary_iterator(u, 1,               2, block_height - 1, 0); break;
    case LEFT:  begin = boundary_iterator(u, 1,               0, 1,                2); break;
    case RIGHT: begin = boundary_iterator(u, block_width - 1, 0, 1,                2); break;
    default:
      CkPrintf("noop send\n");
      return false;
    }
    end = begin + count;

    int k = 0;
    for (; begin != end; ++k, ++begin) {
      boundary[k] = downSample(u, begin.x, begin.y);
      VB(logFile << boundary[k] << "\t";);
    }
    VB(logFile << std::endl;);
    CkAssert(k == count);

    if (which)
      thisProxy(receiver).receiveRefGhosts(iterations, sender_direction, count, boundary);
    else
      thisProxy(receiver).receiveGhosts   (iterations, sender_direction, count, boundary, thisIndex, rand());

    lastSent = receiver;
    return true;
  }

  VB(logFile << thisIndex.getIndexString() << " Will Wait For Ghost from Dir " << dir << ", iteration " << iterations << std::endl;);
  return false;
}

void Advection::begin_iteration(void) {
  //logFile << "String: " << thisIndex.getIndexString() << std::endl;
    
  char fname[100];
  sprintf(fname, "out/out_%s_%d", thisIndex.getIndexString().c_str(), iterations);
  VB(outFile.open(fname););
  VB(logFile << "************************Begin Iteration " << iterations << " on " << thisIndex.getIndexString() << std::endl;);

    for(int i=0; i<3*NUM_NEIGHBORS; i++)
      nbr_dataSent[i]=false;
    
  hasReceived.clear();

  for(int j=1; j<=block_height; j++){
    left_edge[j-1] = u[index(1,j)];
    right_edge[j-1] = u[index(block_width,j)];
  }
  //logFile << "Seding Ghosts " << thisIndex.getIndexString() << std::endl;
  for(int i=0; i<NUM_NEIGHBORS; i++){
    sendGhost(i);
  }
  VB(logFile << "Done Sending Ghosts " << thisIndex.getIndexString() << std::endl;);;
}
template<class T>
void Advection::print_Array(T* array, int size, int row_size){
#ifdef LOGGER
  for(int i=0; i<size; i++){
    if(i%row_size==0)
      logFile << std::endl;
    logFile << array[i] << '\t';
  }
#endif
}

void Advection::process(int iteration, int dir, int size, double gh[]){
//printf("[%d] process %d %d\n", thisIndex, iter, dir);
  VB(logFile << thisIndex.getIndexString() << " received data for direction " << dir << ", iteration " << iteration << ", " << iterations << std::endl;);
    for(int i=0; i<size; i++){
      VB(logFile << gh[i] << '\t';);
        }
    VB(logFile << std::endl;);

    hasReceived.insert(dir);

    boundary_iterator iter;

    if (dir <= RIGHT)
      imsg += 1;
    else
      imsg += 0.5;

    switch(dir){
    case UP:         iter = boundary_iterator(u, 1,                 1, 0,                  0); break;
    case UP_LEFT:    iter = boundary_iterator(u, 1,                 1, 0,                  0); break;
    case UP_RIGHT:   iter = boundary_iterator(u, block_width/2 + 1, 1, 0,                  0); break;
    case DOWN:       iter = boundary_iterator(u, 1,                 1, block_height   + 1, 0); break;
    case DOWN_LEFT:  iter = boundary_iterator(u, 1,                 1, block_height   + 1, 0); break;
    case DOWN_RIGHT: iter = boundary_iterator(u, block_width/2 + 1, 1, block_height   + 1, 0); break;
    case LEFT:       iter = boundary_iterator(u, 0,                 0, 1,                  1); break;
    case LEFT_UP:    iter = boundary_iterator(u, 0,                 0, 1,                  1); break;
    case LEFT_DOWN:  iter = boundary_iterator(u, 0,                 0, block_height/2 + 1, 1); break;
    case RIGHT:      iter = boundary_iterator(u, block_width   + 1, 0, 1,                  1); break;
    case RIGHT_UP:   iter = boundary_iterator(u, block_width   + 1, 0, 1,                  1); break;
    case RIGHT_DOWN: iter = boundary_iterator(u, block_width   + 1, 0, block_height/2 + 1, 1); break;
    default:
      CkAbort("ERROR\n");
    }

    for(int i=0; i<size; ++i, ++iter)
      *iter = gh[i];
}

int simpleDirectionFromComplex(int dir) {
  switch (dir) {
  case UP_RIGHT:   case UP_LEFT:   return UP;
  case DOWN_RIGHT: case DOWN_LEFT: return DOWN;
  case LEFT_DOWN:  case LEFT_UP:   return LEFT;
  case RIGHT_DOWN: case RIGHT_UP:  return RIGHT;
  default: CkAbort("called on non-complex direction");
  }
}

int cornerDirection(int dir) {
  switch (dir) {
  case RIGHT_UP:   return UP_RIGHT;
  case RIGHT_DOWN: return DOWN_RIGHT;
  case LEFT_UP:    return UP_LEFT;
  case LEFT_DOWN:  return DOWN_LEFT;
  case UP_LEFT:    return LEFT_UP;
  case UP_RIGHT:   return RIGHT_UP;
  case DOWN_LEFT:  return LEFT_DOWN;
  case DOWN_RIGHT: return RIGHT_DOWN;
  default: CkAbort("called on other direction");
  }
}

bool Advection::hasReceivedFromDir(int dir) {
  return hasReceived.find(dir) != hasReceived.end();
}

bool Advection::hasReceivedFromAroundCorner(int aroundCorner) {
  return hasReceivedFromDir(aroundCorner) ||
    hasReceivedFromDir(simpleDirectionFromComplex(aroundCorner));
}

void Advection::sendReadyData2RefiningNeighbors(){
  for(int i=NUM_NEIGHBORS; i<3*NUM_NEIGHBORS; i++) {
    if(nbr_decision[i] == REFINE && !nbr_dataSent[i]) {
      if (hasReceivedFromDir(i) &&
          hasReceivedFromAroundCorner(cornerDirection(i))) {
        interpolateAndSend(i);
        nbr_dataSent[i]=true;
      }
    }
  }
}

void Advection::sendReadyData(){
  //check if data can be sent to any of the refined neighbors
  //If the neighbors are at the same level or do not exist at all 
  //data will be sent in begin_iteration function and need 
  //not be sent here
  for(int i=0; i<NUM_NEIGHBORS; i++)
    if(nbr_isRefined[i] && !nbr_dataSent[i]) {
      int sidea, sideb, cornera, cornerb;
      switch(i) {
      case LEFT:  sidea   = LEFT_UP;    sideb   = LEFT_DOWN;
                  cornera = UP_LEFT;    cornerb = DOWN_LEFT;  break;
      case RIGHT: sidea   = RIGHT_UP;   sideb   = RIGHT_DOWN;
                  cornera = UP_RIGHT;   cornerb = DOWN_RIGHT; break;
      case UP:    sidea   = UP_RIGHT;   sideb   = UP_LEFT;
                  cornera = RIGHT_UP;   cornerb = LEFT_UP;    break;
      case DOWN:  sidea   = DOWN_RIGHT; sideb   = DOWN_LEFT;
                  cornera = RIGHT_DOWN; cornerb = LEFT_DOWN;  break;
      }

      if ( hasReceivedFromDir(sidea) && hasReceivedFromDir(sideb) &&
           hasReceivedFromAroundCorner(cornera) && hasReceivedFromAroundCorner(cornerb)
           ) {
        interpolateAndSend(i);
        nbr_dataSent[i]=true;
      }
    }
}

int Advection::getSourceDirection(int NBR) {
  switch (NBR) {
  case UP:    case UP_RIGHT:   case UP_LEFT:    return DOWN;
  case DOWN:  case DOWN_RIGHT: case DOWN_LEFT:  return UP;
  case LEFT:  case LEFT_UP:    case LEFT_DOWN:  return RIGHT;
  case RIGHT: case RIGHT_UP:   case RIGHT_DOWN: return LEFT;
  default: CkAbort("Trying to send in an unknown direction");
  }
}
    
void Advection::interpolateAndSend(int NBR) {
  if(NBR==UP){
    interpolateAndSend(UP_LEFT);
    interpolateAndSend(UP_RIGHT);
    return;
  }
  else if(NBR==DOWN){
    interpolateAndSend(DOWN_LEFT);
    interpolateAndSend(DOWN_RIGHT);
    return;
  }
  else if(NBR==LEFT){
    interpolateAndSend(LEFT_UP);
    interpolateAndSend(LEFT_DOWN);
    return;
  }
  else if(NBR==RIGHT){
    interpolateAndSend(RIGHT_UP);
    interpolateAndSend(RIGHT_DOWN);
    return;
  }

  double *boundary, *out;
  boundary_iterator in, end;
  int count = getGhostCount(NBR);
  int p = 1, m = -1;
  int x1, x2, y1, y2;
  int cdir; // child direction relative to this chare's neighboring uncle
  double sx_r, sx_l, sy_u, sy_d;
  double *sx1, *sx2, *sy1, *sy2;
  switch(NBR) {
  case UP_LEFT:    in = boundary_iterator(u, 1,             1, 1,              0); sx1 = &sx_l; 	sx2 = &sx_r;    sy1 = &sy_u; 	sy2 = &sy_u;  cdir = DOWN_LEFT;  break;
  case UP_RIGHT:   in = boundary_iterator(u, block_width/2, 1, 1,              0); sx1 = &sx_l; 	sx2 = &sx_r;    sy1 = &sy_u; 	sy2 = &sy_u;  cdir = DOWN_RIGHT; break;
  case DOWN_LEFT:  in = boundary_iterator(u, 1,             1, block_height,   0); sx1 = &sx_l; 	sx2 = &sx_r;    sy1 = &sy_d;  sy2 = &sy_d;  cdir = UP_LEFT;    break;
  case DOWN_RIGHT: in = boundary_iterator(u, block_width/2, 1, block_height,   0); sx1 = &sx_l; 	sx2 = &sx_r;    sy1 = &sy_d;  sy2 = &sy_d;  cdir = UP_RIGHT;   break;
  case LEFT_UP:    in = boundary_iterator(u, 1,             0, 1,              1); sx1 = &sx_l; 	sx2 = &sx_l; 	  sy1 = &sy_u; 	sy2 = &sy_d;  cdir = RIGHT_UP;   break;
  case LEFT_DOWN:  in = boundary_iterator(u, 1,             0, block_height/2, 1); sx1 = &sx_l; 	sx2 = &sx_l;    sy1 = &sy_u; 	sy2 = &sy_d;  cdir = RIGHT_DOWN; break;
  case RIGHT_UP:   in = boundary_iterator(u, block_width,   0, 1,              1); sx1 = &sx_r; 	sx2 = &sx_r;    sy1 = &sy_u; 	sy2 = &sy_d;  cdir = LEFT_UP;    break;
  case RIGHT_DOWN: in = boundary_iterator(u, block_width,   0, block_height/2, 1); sx1 = &sx_r; 	sx2 = &sx_r;    sy1 = &sy_u; 	sy2 = &sy_d;  cdir = LEFT_DOWN;  break;
  default: CkAbort("Trying to send to an unrefined or unknown neighbor");
  }
  out = boundary = getGhostBuffer(NBR);
  end = in + count/2;

  for (; in != end; ++in) {
    sx_r = (in.right() - *in) / 4;
    sx_l = -1*(*in - in.left())/4;
    // Possible inversion in definitions of up/down
    sy_u = -1*(*in - in.down())  / 4;
    sy_d = (in.up() - *in)  / 4;
	  
    *out = *in + *sx1 + *sy1; out++;
    *out = *in + *sx2 + *sy2; out++;
  }

  QuadIndex receiver = nbr[simpleDirectionFromComplex(NBR)].getChild(map_child(cdir));
  thisProxy(receiver).receiveGhosts(iterations, getSourceDirection(NBR), count, boundary, thisIndex, rand());

#ifdef LOGGER
	logFile << "sending interpolated data to " << receiver.getIndexString().c_str() << std::endl;
  for(int i=0; i<count; i++){
    logFile << boundary[i] << '\t';
  }
  logFile << std::endl;
	
#endif

  return;
}

void Advection::compute_and_iterate(){
  //logFile << "dt: " << dt << " ap:" << ap << " an:" << an << std::endl;
  //logFile << iterations << std::endl;
  logFile << "entire u before updating" << std::endl;
  for(int j=0; j<block_height+2; j++){
      for(int i=0; i<block_width+2; i++)
          logFile << u[index(i,j)] << "\t";
        logFile << std::endl;
  }

  for(int i=1; i<=block_width; i++){
    for(int j=1; j<=block_height; j++){
      up = (u[index(i+1,j)] - u[index(i,j)])/dx;
      un = (u[index(i,j)]-u[index(i-1,j)])/dx;

      u2[index(i,j)] = u[index(i,j)] - dt* (ap*un + an*up);
      //logFile << u[index(i,j+1)] << ", " << u[index(i,j)] << ", " << u[index(i,j-1)] << std::endl;
      //logFile << thisIndex.getIndexString() << " u2: " << u2[index(i,j)] << std::endl;
    }
  }

  for(int i=1; i<=block_width; i++)
    for(int j=1; j<=block_height; j++){
      up = (u[index(i,j+1)] - u[index(i,j)])/dy;
      un = (u[index(i,j)] - u[index(i,j-1)])/dy;

      u3[index(i,j)] = u[index(i,j)] - dt*(ap*un + an*up);
      //logFile << u[index(i,j+1)] << ", " << u[index(i,j)] << ", " << u[index(i,j-1)] << std::endl;
      //logFile << thisIndex.getIndexString() << "dx " << dx << ", dy " << dy << ", up: " << up << " un: " << un << " u3: " << u3[index(i,j)] << std::endl;
    }

  for(int j=1; j<=block_height; j++)
    for(int i=1; i<=block_width; i++){
      u[index(i,j)] = 0.5*(u2[index(i,j)] + u3[index(i,j)]);
      if(u[index(i,j)] <=1);

    }
    
#ifdef LOGGER
  logFile << "Values of " << thisIndex.getIndexString() << ", iteration " << iterations << std::endl;
  for(int i=1; i<=block_height; i++){
    for(int j=1; j<=block_width; j++)
      logFile << u[index(j,i)] << "\t";
    logFile << std::endl;
  }
  logFile << std::endl;
    
  //outFile << "coordinates: " << xc << ", " << yc << std::endl;
  //outFile << dx << ", " << dy << std::endl;
  for(int i=1; i<=block_width; i++){
    for(int j=1; j<=block_height; j++){
      //outFile << xmin << ", " << double(xc*block_width + i) << std::endl;
      outFile << xmin + (double(i))*dx - 0.5*dx << " "\
              << ymin + (double(j))*dy - 0.5*dy << " "\
              << u[index(i,block_height+1-j)] << std::endl;
    }
  }
  outFile.flush();
  outFile.close();
#endif
}

void Advection::done(){
  if(!has_terminated){
    has_terminated=true;
    //contribute();
    ckout << thisIndex.getIndexString().c_str()  << " is now terminating" << endl;
    if(thisIndex.getDepth()!=min_depth)
      thisProxy(thisIndex.getParent()).done();
  }
}

void Advection::iterate() {
  if(iterations >= max_iterations){
    //ckout << thisIndex.getIndexString().c_str() << " now terminating" << endl;
    VB(logFile << thisIndex.getIndexString() << " now terminating" << std::endl;);
    contribute(CkCallback(CkIndex_Main::terminate(), mainProxy));
    //contribute();
    //if(thisIndex.getDepth()!=min_depth)
    //  thisProxy(thisIndex.getParent()).done();
    return;
  }

  myt = myt+mydt;
  if(myt < tmax){
    mydt = min(dx,dy)/v * cfl;
    if ((myt + mydt) >= tmax )
      mydt = tmax - myt;

    //time to check need for refinement/coarsening
    if(iterations % refine_frequency == 0) {
      VB(logFile << "Entering Mesh Restructure Phase on " << thisIndex.getIndexString() << ", iteration " << iterations << std::endl;);
      //contribute(CkCallback(CkIndex_Advection::startRemesh(), thisProxy));
      startRemesh();
    }
    else {
      VB(logFile << "calling doStep now, iteration " << iterations << std::endl;);
      doStep();
    }
  }
  else {
    CkPrintf("Contribute from %s\n", thisIndex.getIndexString().c_str());

    //contribute();
  }
}

DECISION Advection::getGranularityDecision(){
    delx = 0.5/dx;
    dely = 0.5/dy;
    dely_f = dely;
    double error=0;

    for(int i=1; i <= block_width; i++){
        for(int j=1; j<=block_height; j++){
            // d/dx
            delu[0][i][j] = u[index(i+1, j)] - u[index(i-1,j)];
            delu[0][i][j] = delu[0][i][j]*delx;

            delua[0][i][j] = abs(u[index(i+1, j)]) + abs(u[index(i-1, j)]);
            delua[0][i][j] = delua[0][i][j]*delx;

            // d/dy
            delu[1][i][j] = u[index(i, j+1)] - u[index(i, j-1)];
            delu[1][i][j] = delu[1][i][j]*dely_f;

            delua[1][i][j] = abs(u[index(i, j+1)]) + abs(u[index(i, j-1)]);
            delua[1][i][j] = delua[1][i][j]*dely_f;
        }
    }
    
    int istart=2, iend=block_width-1, jstart=2, jend=block_height-1;
    for (int i=istart;i<=iend;i++){
        for (int j=jstart;j<=jend;j++){

            delu2[0] = (delu[0][i+1][j] - delu[0][i-1][j])*delx;
            delu3[0] = (abs(delu[0][i+1][j]) + abs(delu[0][i-1][j]))*delx;
            delu4[0] = (delua[0][i+1][j] + delua[0][i-1][j])*delx;

            // d/dydx
            delu2[1] = (delu[0][i][j+1] - delu[0][i][j-1])*dely_f;
            delu3[1] = (abs(delu[0][i][j+1]) + abs(delu[0][i][j-1]))*dely_f;
            delu4[1] = (delua[0][i][j+1] + delua[0][i][j-1])*dely_f;

            // d/dxdy
            delu2[2] = (delu[1][i+1][j] - delu[1][i-1][j])*delx;
            delu3[2] = (abs(delu[1][i+1][j]) + abs(delu[1][i-1][j]))*delx;
            delu4[2] = (delua[1][i+1][j] + delua[1][i-1][j])*delx;

            // d/dydy
            delu2[3] = (delu[1][i][j+1] - delu[1][i][j-1])*dely_f;
            delu3[3] = (abs(delu[1][i][j+1]) + abs(delu[1][i][j-1]))*dely_f;
            delu4[3] = (delua[1][i][j+1] + delua[1][i][j-1])*dely_f;
		
            // compute the error
            double num = 0.;
            double denom = 0.;

            for (int kk = 0; kk < ndim2; kk++){  // kk= 1, 2, 3, 4
                num = num + pow(delu2[kk],2.);
                denom = denom + pow(delu3[kk], 2.) + (refine_filter*delu4[kk])*2;
                // refine_filter is the epsilon in my writing sheet
                // refine_filter = 0.01 by default
            }
            // compare the square of the error
            if (denom == 0. && num != 0.){
                error = std::numeric_limits<double>::max();
            }else if (denom != 0.0){
                error = std::max(error, num/denom);
            }
        }
    }
    error = sqrt(error);
    if(error < derefine_cutoff && thisIndex.getDepth() > min_depth)
        return DEREFINE;
    else if(error > refine_cutoff && thisIndex.getDepth() < max_depth)
        return REFINE;  
    else
        return STAY;
}

void Advection::resetMeshRestructureData(){
  /*Reset old data*/
  /*Phase1 Resetting*/
  for(int i=0; i<3*NUM_NEIGHBORS; i++)
    nbr_decision[i] = DEREFINE;//by default a neighbor will derefine
                            //unless indicated otherwise through exchange of a
                            //message in phase1

  decision = INV;

  for(int i=0; i<NUM_CHILDREN; i++)
    child_decision[i]=INV;
  VB(logFile << "setting parentHasAlreadyMadeDecision to false" << std::endl;);
    parentHasAlreadyMadeDecision=false;
  hasReceivedParentDecision=false;
  hasCommunicatedSTAY=false;
  hasCommunicatedREFINE=false;

  /*Phase2 resetting*/
  hasAllocatedMemory=false;
}

void Advection::doRemeshing(){
  //CkPrintf("%s doMeshRestructure %d\n", thisIndex.getIndexString().c_str(), iterations);

  if(!hasReset){
    hasReset=true;
    resetMeshRestructureData();
  }
	
  if(!isRefined) {//run this on leaf nodes
    /*Done Resetting Old Data*/
    /* It is possible that some Message has alreasdy been processed
       MeshRestructure Phase was called*/
    VB(logFile << thisIndex.getIndexString() << " decision before getGranularityDecision is " << decision << std::endl;);

    if (decision == REFINE);//no need to call getGranularityDecision
    else
        decision = max(decision, getGranularityDecision());


    //ckout << thisIndex.getIndexString().c_str() << " decision = " << decision << ", iteration = " << iterations << endl;
    //initiate Phase1 of the computation
    updateNeighborsofChangeInDecision();
    //Inform Parent About Start of the Restructure Phase 
    // Think of a better way to do this because right now as many as 
    // number of unrefined children are informing the parent about the
    // start of the restructure phase while only one of them need to do so
    if(parent!=thisIndex) {
      //CkPrintf("%s -> %s startRemesh %d\n", thisIndex.getIndexString().c_str(),
      //parent.getIndexString().c_str(), iterations);
      //thisProxy(parent).doMeshRestructure();
    }
  }
  else if(isGrandParent && !parentHasAlreadyMadeDecision){
    //CkPrintf("%s was a grandparent %d\n", thisIndex.getIndexString().c_str(), iterations);
    informParent(-1, INV);
  } else {
    //CkPrintf("%s was a parent %d\n", thisIndex.getIndexString().c_str(), iterations);
    //thisProxy[thisIndex].startRemesh();
  }
}

/***** PHASE1 FUNCTIONS****/
void Advection::communicatePhase1Msgs(){
    
  if(decision == REFINE || decision == STAY){
    //tell the neighbor if it exists, also tell to children of the neighbor 
    // if that neighbor is refined
    //In case the neighbor does not exist, just tell to the parent of the 
    //non-existing neighbor

    for(int i=0; i<NUM_NEIGHBORS; i++){
      if(nbr_exists[i] && !nbr_isRefined[i]){
        VB(logFile << thisIndex.getIndexString() << " sending decision " << decision << " to " << nbr[i].getIndexString() << std::endl;);
          thisProxy(nbr[i]).exchangePhase1Msg(SENDER_DIR[i], decision);//Since Phase1Msgs are only refinement messages
      }
      //just send your direction w.r.t. to the receiving neighbor
      else if(nbr_exists[i] && nbr_isRefined[i]){
        //Get Corresponding Children of the neighbor
        QuadIndex q1, q2;
        getChildren(nbr[i], SENDER_DIR[i], q1, q2);
        VB(logFile << thisIndex.getIndexString() << " sending decision to " << q1.getIndexString() << std::endl;);
        VB(logFile << thisIndex.getIndexString() << " sending decision to " << q2.getIndexString() << std::endl;);
                    
                    
        thisProxy(q1).exchangePhase1Msg(SENDER_DIR[i], decision);
        thisProxy(q2).exchangePhase1Msg(SENDER_DIR[i], decision);
      }
      else{//send to the parent of the non-existing neighbor
        VB(logFile << thisIndex.getIndexString() << " sending decision " << decision << " to " << nbr[i].getParent().getIndexString() << std::endl;);
        thisProxy(nbr[i].getParent()).exchangePhase1Msg(map_nbr(thisIndex.getQuadI(), i), decision);
      }
    }
  }
  else{//No Need to do anything, just wait for 
    //neighbors to tell if they wish to derefine
  }

  //If my DECISION is to stay or to REFINE, tell the parent
  if((decision == REFINE || decision == STAY) && (parent != thisIndex)){
    thisProxy(parent).informParent(thisIndex.getChildNum(), decision);
  }
}

void Advection::informParent(int childNum, DECISION dec){//Will be called from two contexts: 
  //a) If the parent is a grandparent and also have children that are leaves
  //b) when a child sends REFINE/STAY message to the parent
  if(!hasReset){
    hasReset=true;
    resetMeshRestructureData();
  }
  VB(logFile << thisIndex.getIndexString() << ": in informParent called by child " << childNum << " with decision = " << dec << ", iteration " << iterations << std::endl;);
  if(childNum >= 0)
      child_decision[childNum]=dec;

  if(dec==REFINE){
    child_isRefined[childNum]=true;
    isGrandParent = true;
  }
  if(parentHasAlreadyMadeDecision == false){
    VB(logFile << "settin parentHasAlreadyMadeDecision to true, iterations " << iterations << std::endl;);
      parentHasAlreadyMadeDecision = true;
    //tell rest of the children which are not refined
    for(int i=0 ;i<NUM_CHILDREN; i++){
      if(i!=childNum && !child_isRefined[i]) {
        thisProxy(thisIndex.getChild(i)).recvParentDecision();
      }
    }
    //inform your parent that you are not going to derefine
    if(parent!=thisIndex)
        thisProxy(parent).informParent(thisIndex.getChildNum(), STAY);
  }
}

void Advection::recvParentDecision(){
  VB(logFile << thisIndex.getIndexString() << " has received decision from parent " << std::endl;);
    if(!hasReset){
      hasReset=true;
      resetMeshRestructureData();
    }

  hasReceivedParentDecision = true;
  decision = std::max(STAY, decision);
  if(!isRefined)
    updateNeighborsofChangeInDecision();
}

bool isDirectionSimple(int dir) {
  return dir == LEFT || dir == RIGHT || dir == UP || dir == DOWN;
}

void Advection::updateNeighborsofChangeInDecision(){
  if(decision == REFINE && !hasCommunicatedREFINE){
    hasCommunicatedREFINE=true;
    communicatePhase1Msgs();
  }else if(decision == STAY && !hasCommunicatedSTAY){
    hasCommunicatedSTAY=true;
    communicatePhase1Msgs();
  }
}

void Advection::exchangePhase1Msg(int dir, DECISION remoteDecision){//Phase1 Msgs are either REFINE or STAY messages
  VB(CkAssert((remoteDecision == REFINE || remoteDecision == STAY)););
  VB(logFile << thisIndex.getIndexString() << " received decision " << remoteDecision << " from direction " << dir << std::endl; );
  if(!hasReset){
    hasReset=true;
    resetMeshRestructureData();
  }
  
  nbr_decision[dir] = std::max(remoteDecision, nbr_decision[dir]);
  remoteDecision = nbr_decision[dir];

    /*inv=-1, coarsen=0, stay=1, refine=2
    refine/stay messages with remoteDecision
    if sender is sibling  //note that this situation will be handled by the informParent and recvParentDecision methods
        myDecision = max(myDecision, remoteDecision)
    if sender is from simple direction
        if remoteDecision == stay
            myDecision = myDecision//I can do whatever I want to do
        else if remoteDecision == refine
            myDecision = max(stay, myDecision)
    else if sender is uncle
        myDecision = myDecision // I can do whatever I want to do
    else if sender is a refined neighbor
        myDecision = max(myDecision, remoteDecision)*/
  //logFile << nbr_exists[dir] << ", " << nbr_isRefined[dir] << std::endl;  
  if(isDirectionSimple(dir) && nbr_exists[dir]){
    //logFile << "i am here1, remoteDecision = " << remoteDecision << ", " << decision << std::endl;
    if(remoteDecision == STAY);
        //decision=std::max(STAY, decision);
    else if(remoteDecision == REFINE)
        decision = std::max(STAY, decision);
    //decision = std::max(STAY, decision);
  }
  else if (isDirectionSimple(dir) && !nbr_exists[dir])
  /*{
      logFile << "i am here2" << std::endl;
        decision=max(STAY, decision);*/
    ;//decision = std::max(STAY, decision);
  /*}*/
  else if(!isDirectionSimple(dir)){
    //logFile << "i am here3" << std::endl;
    VB(CkAssert(isDirectionSimple(dir) == false););
    decision = std::max(decision, remoteDecision);      

  }
  else{
      CkAbort("unacceptable condition");
  }


  VB(logFile << thisIndex.getIndexString() << " decision: " << decision << std::endl;);
  updateNeighborsofChangeInDecision();
}

#define index_c(i,j) (int)((j)*(block_width/2) + i)
ChildDataMsg::ChildDataMsg(bool isInMeshGenerationPhase, int cnum, double mt, double mdt, int meshGenIterations, int iter, double* u, bool* nbr_exists, bool* nbr_isRefined, DECISION* nbr_decision){
  //logFile << child_nbr_isRefined << std::endl;
  if(isInMeshGenerationPhase==false){
    for(int i=1; i<= block_width; i+=2){
      for(int j=1; j<=block_height; j+=2){
        int idx = index_c(i/2, j/2);
        child_u[idx] = downSample(u, i, j);
        //logFile << child_u[idx] << "\t";
      }
      //logFile << std::endl;
    }
  }
  this->isInMeshGenerationPhase = isInMeshGenerationPhase;
  this->meshGenIterations = meshGenIterations;
  childNum = cnum;
  iterations=iter;
  myt=mt;
  mydt=mdt;
  //logFile << child_nbr_exists << ", " << child_nbr_isRefined << ", " << child_nbr_decision << std::endl;
  memcpy(child_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
  //logFile << NUM_NEIGHBORS << std::endl;
  memcpy(child_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
  memcpy(child_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);

}

void getRefinedNbrDirections(int dir, int &d1, int &d2){//returns the direction numbers of the refined neighbors in direction 'dir'
  switch(dir){
    case RIGHT: d1=RIGHT_UP;  d2=RIGHT_DOWN; return;
    case LEFT:  d1=LEFT_UP;   d2=LEFT_DOWN; return;
    case UP:    d1=UP_LEFT;   d2=UP_RIGHT; return;
    case DOWN:  d1=DOWN_LEFT; d2=DOWN_RIGHT; return;
  }
}

/**** PHASE2 FUNCTIONS ****/
void Advection::doPhase2(){
  //cout << thisIndex.getIndexString() << " starting phase 2 " << iterations << std::endl;
  VB(logFile << thisIndex.getIndexString() << " Entering Phase2, iteration " << iterations << std::endl;);

  VB(logFile << thisIndex.getIndexString() << " decision = " << decision << std::endl;);

  /*
  //Send Appropriate Data to the Neighbors based on their New Status
  if(isRefined==false) {// can send data only if I am a leaf node
    for(int i=0; i<NUM_NEIGHBORS; i++) {//check if any of the neighbors need data for its refinement
      //if(nbr_decision[i]==REFINE){
      if(((nbr_exists[i] &&!nbr_isRefined[i])||!nbr_exists[i]) && nbr_decision[i]==REFINE){
        if(i==LEFT){
          for(int j=1; j<=block_height; j++)
            left_edge[j-1] = u[index(1,j)];
        }
        else if(i==RIGHT){
          for(int j=1; j<=block_height; j++)
            right_edge[j-1] = u[index(block_width,j)];
        }

        //to serve getGhostsAndRefine()
        if(sendGhost(i))
          ;
      }else{
        //wait for data to arrive and then send it back
        int d1, d2;
        if(i==RIGHT){d1=RIGHT_UP; d2=RIGHT_DOWN;}
        else if(i==LEFT){d1=LEFT_UP; d2=LEFT_DOWN;}
        else if(i==UP){d1=UP_LEFT; d2=UP_RIGHT;}
        else if(i==DOWN){d1=DOWN_LEFT; d2=DOWN_RIGHT;}
        if(nbr_decision[d1]==REFINE) {
          thisProxy[thisIndex].getAndSendGhost();
        }
        if(nbr_decision[d2]==REFINE) {
          thisProxy[thisIndex].getAndSendGhost();
        }
      }
    }

    //If I am going to refine check if any of the neighbor needs my ghosts for it to send its boundaries to me
    if(decision==REFINE)
      for(int i=0; i<NUM_NEIGHBORS; i++){
        if(!nbr_exists[i] && nbr_decision[i]!=REFINE){//send the ghost layer; dont send if the neighbors decision is refine because in that case we
          //already sent the ghost layer in the loop above
          if(i==LEFT){
            for(int j=1; j<=block_height; j++)
              left_edge[j-1] = u[index(1,j)];
          }
          else if(i==RIGHT){
            for(int j=1; j<=block_height; j++)
              right_edge[j-1] = u[index(block_width,j)];
          }
          //to serve getGhostsAndRefine()
          if(sendGhost(i))
            ;
        }
      }
  }*/

  //CkPrintf("%s middle of phase 2 %d\n", thisIndex.getIndexString().c_str(), iterations);

  if(isRefined && !isGrandParent && !parentHasAlreadyMadeDecision){//I am a parent(whose None of the Children Are Refined) and has to derefine
    //Get Data From the Children and extrapolate it
  }
  else if(decision == DEREFINE && !isRefined){//send data to the parent
    VB(logFile << thisIndex.getIndexString() << " Sending Values to Parent" << std::endl;;);
    size_t sz = ((block_height)*(block_width))/4;
    ChildDataMsg *msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS) 
                            ChildDataMsg(0, thisIndex.getChildNum(), myt, mydt, meshGenIterations, iterations, u, nbr_exists, nbr_isRefined, nbr_decision);

    thisProxy(parent).recvChildData(msg);
    //deallocate all your memory and destroy yourself
    VB(logFile << "Destroying " << thisIndex.getIndexString() << std::endl;);
    shouldDestroy = true;
    thisProxy[thisIndex].ckDestroy();
    VB(logFile << "Done Destroying " << thisIndex.getIndexString() << std::endl;);
  }
  else if(decision == REFINE){
    VB(logFile << "Refine called on " << thisIndex.getIndexString() << std::endl;);
    //thisProxy[thisIndex].getGhostsAndRefine();
    refine();
  }
   
  updateNbrStatus();

  //logFile << thisIndex.getIndexString() << " decision = " << decision << std::endl;
  VB(logFile << "setting parentHasAlreadyMadeDecision to false" << endl;);
  parentHasAlreadyMadeDecision = false;
  hasReset=false;

}

void Advection::updateNbrStatus(){
  //Update the Status of Your Neighbors, need to be done only if you are going to stay in that position
  if(!isRefined && decision == STAY){
    VB(logFile << "Phase2: " << thisIndex.getIndexString() << " updating the Status of Neighbors" << std::endl;);
      for(int i=0; i<NUM_NEIGHBORS; i++){
        if(!nbr_exists[i]){
          switch(nbr_decision[i]){
            case DEREFINE:
              logFile << "ERROR(" << thisIndex.getIndexString() << "): doPhase(): Uncle(" << i << ") Cannot DEREFINE while I want to STAY" << std::endl;
              ckout << thisIndex.getIndexString().c_str() << " 's nbr in direction " << i << " is set to derefine" << endl; 
              CkAbort("undefined state");
            case REFINE: nbr_exists[i]=true; nbr_isRefined[i]=false; break;
            case STAY: nbr_exists[i]=false; break;
            default: CkAbort("nbr_decision not set");
          }
        }
        else if(nbr_exists[i] && !nbr_isRefined[i]){
          switch(nbr_decision[i]){
            case REFINE: nbr_isRefined[i]=true; break;
            case STAY: nbr_exists[i]=true; nbr_isRefined[i]=false; break;
            case DEREFINE: nbr_exists[i]=false; break;
            default: CkAbort("nbr_decision not set");
          }
        }
        else if(nbr_exists[i] && nbr_isRefined[i]){
          int d1, d2;
          getRefinedNbrDirections(i, d1, d2);
          
          switch(nbr_decision[d1]){
            case STAY:    nbr_exists[i]=true; nbr_isRefined[i]=true;  break;
            case DEREFINE: nbr_exists[i]=true; nbr_isRefined[i]=false; break;
            case REFINE:  CkAbort("unacceptable decision");
            default:      CkAbort("nbr_decision not set");
          }
        }
      }
  }
  else if(decision == REFINE){//I will Now become Inactive and therefore I need not Store Neighbor Status
    isRefined=true;
    isGrandParent=false;
  }else if(isRefined && !isGrandParent && !parentHasAlreadyMadeDecision){// parent going to destroy its children
    isRefined = false;
  }
  if(isGrandParent){
    for(int i=0; i<NUM_CHILDREN; i++){
        if(child_decision[i]==INV && child_isRefined[i])//did not receive any message
            child_isRefined[i]=false;
    }
    isGrandParent=false;
    for(int i=0; i<NUM_CHILDREN; i++)
        if(child_isRefined[i])
            isGrandParent=true;
  }
  VB(logFile << "isGrandparent = " << isGrandParent << ", iteration = " << iterations << std::endl;);

}

void Advection::recvChildData(ChildDataMsg *msg){
  VB(logFile << "Mem Check at Beginning of recvChildData" << std::endl;);
  VB(logFile << thisIndex.getIndexString() << " received data from Child " << msg->childNum << " for coarsening" << std::endl;);
    myt = msg->myt;
  mydt = msg->mydt;
  iterations = msg->iterations;
  meshGenIterations = msg->meshGenIterations;

  if(!hasAllocatedMemory){
    hasAllocatedMemory=true;
    mem_allocate_all();
  }
  int st_i, end_i, st_j, end_j;
  switch(msg->childNum){
      case 0:  st_i = block_width/2+1;  end_i = block_width;    st_j = 1;                 end_j = block_height/2;  break;
      case 1:  st_i = 1;                end_i = block_width/2;  st_j = 1;                 end_j = block_height/2;  break;
      case 2:  st_i = 1;                end_i = block_width/2;  st_j = block_height/2+1;  end_j = block_height;    break;
      case 3:  st_i = block_width/2+1;  end_i = block_width;    st_j = block_height/2+1;  end_j = block_height;    break;
      default: CkAbort("undefined child number");
  }
    
  int ctr=0; double rsq;
  if(msg->isInMeshGenerationPhase)
      applyInitialCondition();
  else{
    for(int j=st_j; j<=end_j; j++){
      for(int i=st_i; i<=end_i; i++){
        u[index(i,j)]=msg->child_u[ctr];
        VB(logFile << msg->child_u[ctr] << ", " << u[index(i,j)] << "\t";);
        ctr++;
      }
      VB(logFile << std::endl;);
    }
  }

  //Update the Status of Your Neighbors based on Data Sent from the Children
  int c1, c2;
  switch(msg->childNum){
    case 0: c1=RIGHT; c2=UP;   break;
    case 1: c1=LEFT;  c2=UP;   break;
    case 2: c1=LEFT;  c2=DOWN; break;
    case 3: c1=RIGHT; c2=DOWN; break;
  }
    
  this->iterations=msg->iterations;
  this->myt=msg->myt;
  this->mydt=msg->mydt;
  setNbrStatus(c1, msg);
  setNbrStatus(c2, msg);
  
  decision=DEREFINE;
  delete msg;
}

inline void Advection::setNbrStatus(int dir, ChildDataMsg* msg){
  if(msg->child_nbr_exists[dir]){
    nbr_exists[dir]=true;
    if(msg->child_nbr_isRefined[dir])
      nbr_isRefined[dir]=true;
    if(!msg->child_nbr_isRefined[dir] && msg->child_nbr_decision[dir]==DEREFINE){
      nbr_isRefined[dir]=false;
    }
    if(!msg->child_nbr_isRefined[dir] && msg->child_nbr_decision[dir]==STAY){
      nbr_isRefined[dir]=true;
    }
  }
  else if(!msg->child_nbr_exists[dir]){
    if(msg->child_nbr_decision[dir]==REFINE){
      nbr_exists[dir]=true;
      nbr_isRefined[dir]=true;
    }
    else if(msg->child_nbr_decision[dir]==STAY){
      nbr_exists[dir]=true;
      nbr_isRefined[dir]=false;
    }
    else {//DEREFNE
      nbr_exists[dir]=false;
    }
  }    
}

void Advection::requestNextFrame(liveVizRequestMsg *m){
  if (isRefined){
    liveVizDeposit(m, 0, 0, 0, 0, NULL, this);//submit a zero size array
  }

  int sy = yc;//thisIndex%num_chare_cols;
  int sx = xc;//thisIndex/num_chare_rows;
  int w= block_width, h = block_height;

  unsigned char *intensity = new unsigned char[3*w*h];
  for(int j=0; j<h; j++){
    for(int i=0; i<w; i++){
      // logFile << u[i+1][j+1] << std::endl;
      intensity[3*(i+j*w)+ 0] = 255;// red component
      intensity[3*(i+j*w)+ 1] = 255 - (u[index(i+1,j+1)]-1)*255;// BLUE component
      intensity[3*(i+j*w)+ 2] = 255 - (u[index(i+1,j+1)]-1)*255;// green component
    }
  }
  liveVizDeposit(m, sy*w,sx*h, w,h, intensity, this);
  delete[] intensity;
}
#define index_l(i,j)  (int)((j)*(block_width) + i)

void Advection::interpolate(double *u, vector<double>& refined_u, int xstart, int xend, int ystart, int yend){
  double sx_l, sx_r, sy_u, sy_d;
  int m=1, n=1;
  for(int i=xstart; i<=xend; i++){
    for(int j=ystart; j<=yend; j++){
      sx_l = (-u[index(i-1,j)]+u[index(i,j)]) / 4;
      sx_r = (-u[index(i,j)]+u[index(i+1,j)]) / 4;

      sy_u = (-u[index(i,j-1)]+u[index(i,j)]) / 4;
      sy_d = (-u[index(i,j)]+u[index(i,j+1)]) / 4;

      refined_u[index_l(2*(m-1),   2*(n-1))  ] = u[index(i,j)] - sx_l - sy_u;
      refined_u[index_l(2*(m-1)+1, 2*(n-1))  ] = u[index(i,j)] + sx_r - sy_u;
      refined_u[index_l(2*(m-1),   2*(n-1)+1)] = u[index(i,j)] - sx_l + sy_d;
      refined_u[index_l(2*(m-1)+1, 2*(n-1)+1)] = u[index(i,j)] + sx_r + sy_d;
		  VB(logFile << u[index(i,j)] << ", " << sx_l << ", " << sx_r << ", " << sy_u << ", " << sy_d << std::endl;);
      n++;
    }
    m++;
    n=1;
  }
}

void Advection::refineChild(unsigned int sChild, int xstart, int xend, int ystart, int yend, double xmin, double ymin) {
  QuadIndex child = thisIndex.getChild(sChild);

  size_t sz = (block_width)*(block_height);
  vector<double> refined_u(sz);
	
	VB(logFile << "interpolating data for child: " << child.getIndexString().c_str() << std::endl;);
  interpolate(u, refined_u, xstart, xend, ystart, yend);

  InitRefineMsg * msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)
    InitRefineMsg(0, dx/2, dy/2, myt, mydt, xmin, ymin, meshGenIterations, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
  thisProxy(child).insert(msg);
}

void Advection::refine(){
  //Spawn The four children and give them the data
  //Assuming we already have the new boundary data
    
  //Interpolate the data and give it to the children when they are initialized
  // boundaries of the children will have to be sent by the neighbor
  VB(logFile << thisIndex.getIndexString() << " is refining" << std::endl;);

  refineChild(1, 1,               block_width/2, 1,                block_height/2, xmin,           ymin+(ny*dy)/2);
  refineChild(0, block_width/2+1, block_width,   1,                block_height/2, xmin+(nx*dx)/2, ymin+(ny*dy)/2);
  refineChild(2, 1,               block_width/2, block_height/2+1, block_height,   xmin,           ymin);
  refineChild(3, block_width/2+1, block_width,   block_height/2+1, block_height,   xmin+(nx*dx)/2, ymin);

  thisProxy.doneInserting();

  //CkPrintf("%s parent phase 2b iteration %d\n", thisIndex.getIndexString().c_str(), iterations);
  //fflush(stdout);
//   contribute(CkCallback(CkReductionTarget(Advection, phase2Done), thisProxy));

  VB(logFile << thisIndex.getIndexString() << " done with refinement" << std::endl;);;
}

Advection::Advection(InitRefineMsg* msg)
  /*: AdvTerm(thisProxy, thisIndex, true), CBase_Advection()*/
{
  //ckout << thisIndex.getIndexString().c_str() << " created 2" << endl;
  __sdag_init();
  hasReset = false;
  shouldDestroy = false;
  //rootTerminated();
  usesAtSync=CmiTrue;
  has_terminated=false;

  char fname[100];
  sprintf(fname, "log/%s.log", thisIndex.getIndexString().c_str());

  VB(logFile.open(fname););
    //srand(thisIndex.getQuadI() + atoi(thisIndex.getIndexString()));

//Called as a result of refinement of parent
  VB(logFile << "Inserting New Zone: " << thisIndex.getIndexString() << std::endl;);
  this->isRefined = false;

  for(int dir=UP; dir<=RIGHT; ++dir)
    nbr[dir] = thisIndex.getNeighbor(dir);
  if(thisIndex.nbits!=0)
    parent = thisIndex.getParent();

  isGrandParent = false;
  /*set the status of the neighbors
    1. If parent of the neighbor is same as mine, 
    then the status is also the same, 
    i.e. it exists and is not refined 
    because it has just been created alongwith me.
    2. If parent of the neihgbor is not the same as mine
    then infomration about the neighbor will be obtained in two steps
    1. ask the corresponding neighbor of the parent to know if it is refined
    2. if it is refined ask the neighbor if it is refined.
  */
  for(int dir=0; dir<NUM_NEIGHBORS; dir++){
    VB(logFile << thisIndex.getIndexString() << " neighbor in direction " << dir << " is " << nbr[dir].getIndexString() << std::endl;);
      if(nbr[dir].getParent() == thisIndex.getParent()){//if parents are the same, the neighbor has also just been created
                                                        // so it can not be refined
        nbr_exists[dir]=true;
        nbr_isRefined[dir]=false;
      }
      else{//when the parents are not the same
        //if I want to refine then no-one can stop me from doing so
        //so it is possible that my parents neighbor do not exist
        //at this moment but a notification has been sent that they
        //should be generated
        if (msg->parent_nbr_exists[dir] && !msg->parent_nbr_isRefined[dir]){
          switch(msg->parent_nbr_decision[dir]){
            case REFINE: nbr_exists[dir]=true;  nbr_isRefined[dir]=false; break;
            case STAY:   nbr_exists[dir]=false; break;
            case DEREFINE: CkAbort("this neighbor cannot derefine");
            default: VB(logFile << thisIndex.getIndexString() << " nbr decision not set" << endl;); 
                CkAbort("nbr decision not set");
          }
        }
        else if (msg->parent_nbr_exists[dir] && msg->parent_nbr_isRefined[dir]){
            int nbr_dir_wrt_parent = getNbrDir(thisIndex.getQuadI(), dir);//neighbor direction w.r.t. the parent
            switch(msg->parent_nbr_decision[nbr_dir_wrt_parent]){
                case DEREFINE: nbr_exists[dir]=false; break;
                case STAY: nbr_exists[dir]=true; nbr_isRefined[dir]=false; break;
                case REFINE: nbr_exists[dir]=true; nbr_isRefined[dir]=true; break;
                default: CkAbort("nbr decision not set");
            }
        }
        else if (!msg->parent_nbr_exists[dir]){
          VB(CkAssert(msg->parent_nbr_decision[dir]==REFINE););
          nbr_exists[dir]=false;
        }
        else{
            CkAbort("huh.. we should never reach here");
        }
      }
  }

  //Now initialize xmin, xmax, ymin, ymax, dx, dy, myt, mydt
  hasReset=false;
  dx = msg->dx;
  dy = msg->dy;

  myt = msg->myt;
  mydt = msg->mydt;

  xmin = msg->xmin;
  ymin = msg->ymin;
    
  nx = array_height/(num_chare_cols);
  ny = array_width/(num_chare_rows);


  VB(logFile << "xmin: " << xmin << ", ymin: " << ymin << std::endl;);

  thisIndex.getCoordinates(xc, yc);
  meshGenIterations = msg->meshGenIterations;
  iterations = msg->iterations;

  mem_allocate_all();

  for (int i = 0; i < NUM_CHILDREN; ++i)
    child_isRefined[i] = false;
    
  //x[i] and y[i] need not be initialized they are needed only for setting initial conditions
  //delete [] x;
  //delete [] y;
 	
  if(msg->isInMeshGenerationPhase){//setup x[i] and y[i] and initial values
    for(int i=0; i<block_width+2; i++)
      x[i] = xmin + double(i)*dx - 0.5*dx;

    for(int i=0; i<block_height+2; i++)
      y[i] = ymin + double(i)*dy - 0.5*dy;
    
    applyInitialCondition();
  }else{
      //Initialize u - For boundaries I have to wait for the neighbors
      //to send the values, rest of it can be initialized by the values 
      //received from the parent
      int ctr=0;
      for(int j=1; j<=block_height; j++)
        for(int i=1; i<=block_width; i++)
          u[index(i,j)]=msg->refined_u[ctr++];
  }

#ifdef LOGGER
  logFile << "New Child Values" << std::endl;
  for(int i=0; i<block_height; i++){
    for(int j=0; j<block_width; j++)
      logFile << u[index(j+1,i+1)] << "\t";
    logFile << std::endl;
  }
#endif
  decision=REFINE;//to be used for reduction while doing mesh generation
                  // to keep track if any change happened in the last mesh gen iteration 
  //delete the message
  delete msg;

  //call Quiesence detection to begin next iteration*/
  //CkStartQD(CkIndex_Advection::doStep(), &thishandle);
}

void Advection::updateMesh(){
    if(!isRefined){
        //ckout << thisIndex.getIndexString().c_str() << " decision: " << decision << endl;
        if(decision==REFINE){
          double cxmin, cymin;
          QuadIndex child;
          std::vector<double> eVector; //empty vector
          for(int childNum=0; childNum<NUM_CHILDREN; childNum++){
            switch(childNum){
              case 0:  cxmin = xmin + (nx*dx)/2; cymin = ymin + (ny*dy)/2; break;
              case 1:  cxmin = xmin;             cymin = ymin + (ny*dy)/2; break;
              case 2:  cxmin = xmin;             cymin = ymin;             break;
              case 3:  cxmin = xmin + (nx*dx)/2; cymin = ymin;             break;
              default: CkAbort("this is impossible");
            }
            child = thisIndex.getChild(childNum);
            InitRefineMsg * msg = new (0, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)
                                InitRefineMsg(1, dx/2, dy/2, myt, mydt, cxmin, cymin, meshGenIterations, iterations, eVector, nbr_exists, nbr_isRefined, nbr_decision);
            thisProxy(child).insert(msg);
          }
          //ckout << thisIndex.getIndexString().c_str() << " is now refining" << endl;
          thisProxy.doneInserting();
        }
        else if(decision==DEREFINE){
          ChildDataMsg *msg = new (0, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS) 
                                    ChildDataMsg(1, thisIndex.getChildNum(), myt, mydt, meshGenIterations, iterations, 
                                                    u, nbr_exists, nbr_isRefined, nbr_decision);

          thisProxy(parent).recvChildData(msg);
          //deallocate all your memory and destroy yourself
          VB(logFile << "Destroying " << thisIndex.getIndexString() << std::endl;);
          shouldDestroy = true;
          //ckout << thisIndex.getIndexString().c_str() << " is now coarsening" << endl;
          thisProxy[thisIndex].ckDestroy();
        }

    }

    updateNbrStatus(); //update the status of your neighbors
    if(decision==REFINE || decision==DEREFINE)
        decision=INV;
    VB(logFile << "setting parentHasAlreadyMadeDecision to false" << endl;);
    parentHasAlreadyMadeDecision = false;
    hasReset = false;
}

void Advection::ResumeFromSync(){
    phase2Done();
}

void Advection::startLdb(){
    if(thisIndex.bitVector==0)
        ckout << "starting load balancing now.." << endl;
    AtSync();
}
#include "Advection.def.h"
