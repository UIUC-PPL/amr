#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "charm++.h"
#include "liveViz.h"
#include "trace-projections.h"
#include "Main.h"
//#include <algorithm>
extern CProxy_Main mainProxy;

extern int array_height;
extern int array_width;

extern int num_chare_rows;
extern int num_chare_cols;

extern int block_width, block_height;

extern int min_depth, max_depth;

#define wrap_x(a) (((a)+num_chare_rows)%num_chare_rows)
#define wrap_y(a) (((a)+num_chare_cols)%num_chare_cols)


extern CkArrayID a;

extern int nframe;
extern double xctr, yctr, radius;
extern double v;
extern double ap, an;
extern double tmax, t, dt, cfl;
extern map<string, DIR> nbrDirectionMap;
extern map<DIR, DIR> reverse_dir_map;
extern int max_iterations, refine_frequency;

#define index(i,j)  (int)((j)*(block_width+2) + i)

double ***delu, ***delua;
double *delu2, *delu3, *delu4;

PerProcessorChare::PerProcessorChare(){
  delu = new double**[2];
}

InitRefineMsg::InitRefineMsg(double dx, double dy, double myt, double mydt,
                             double xmin, double  ymin, int iterations_,
                             vector<double>& refined_u, bool *nbr_exists,
                             bool *nbr_isRefined, DECISION *nbr_decision) {
  this->dx = dx;
  this->dy = dy;
  this->myt = myt;
  this->mydt = mydt;
  this->xmin = xmin;
  this->ymin = ymin;
  this->iterations = iterations_;
  memcpy(this->refined_u, &refined_u[0], sizeof(double)*block_height*block_width);
  memcpy(this->parent_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
  memcpy(this->parent_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
  memcpy(this->parent_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);
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
  : AdvTerm(thisProxy, thisIndex, true)
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
    VB(logFile << "exists: " << exists << std::endl;
       logFile << "isRefined: " << isRefined << std::endl;
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
    /*for(int i=0; i<block_width+2; i++){
      for(int j=0; j<block_height+2; j++){
      rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[j]-yctr);
      if(i==block_width/2)
      u[index(i, block_height+1-j)] = 1;
      else u[index(i, block_height+1-j)] = 1;
      }
      }*/
    for(int i=0; i<block_width+2; i++){
      for(int j=0; j<block_height+2; j++){
        rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[j]-yctr);
        if(rsq <= radius*radius)
          u[index(i, block_height+1-j)] = 2;
        else u[index(i, block_height+1-j)] = 1;
      }
    }
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
  for(int i=0; i<NUM_NEIGHBORS; i++){
    p|nbr_exists[i];
    p|nbr_isRefined[i];
  }

  for(int i=0; i<3*NUM_NEIGHBORS; i++)
    p|nbr_dataSent[i];

  for(int i=0; i<NUM_CHILDREN; i++)
    p|child_isRefined[i];

  p|parent;
  p|xc;
  p|yc;
  p|imsg;

  if(p.isUnpacking())
    mem_allocate_all();
    
  for(int i=0; i<(block_width+2)*(block_height+2); i++){
    p|u[i];
    p|u2[i];
    p|u3[i];
  }

  for (int i=0; i<block_width+2; i++){
    p|x[i];
    p|top_edge[i];
    p|bottom_edge[i];
  }

  for (int i=0; i<block_height+2; i++){
    p|y[i];
    p|left_edge[i];
    p|right_edge[i];
  }

  p|iterations;
  p|up;
  p|un;
  p|myt;
  p|mydt;
    
  p|hasReceived;
}
    
Advection::~Advection(){
  CkPrintf("%s destructor %d\n", thisIndex.getIndexString().c_str(), iterations);
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

bool Advection::sendGhost(int dir, bool which=0){
  if(nbr_exists[dir] && !nbr_isRefined[dir]){
    int val = rand();
    VB(logFile << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Neighbor " << nbr[dir].getIndexString() << ", iteration " << iterations << ", " << SENDER_DIR[dir] << ", " << val << std::endl;);

    int count;
    double* boundary;
    switch (dir) {
    case LEFT:  count = block_height; boundary = left_edge;                  break;
    case RIGHT: count = block_height; boundary = right_edge;                 break;
    case UP:    count = block_width;  boundary = &u[index(1, 1)];            break;
    case DOWN:  count = block_width;  boundary = &u[index(1, block_height)]; break;
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
    VB(logFile << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Uncle: " << receiver.getIndexString() << ", iteration " << iterations << endl;);
      if(dir==LEFT) {
        for(int j=1; j<=block_height; j+=2){
          left_edge[j/2] = (u[index(1,j)] + u[index(2,j)] + u[index(1,j+1)] +u[index(2,j+1)])/4;
          VB(logFile << left_edge[j/2] << "\t";);
            }
        VB(logFile << std::endl;);
          if(!which)
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getQuadI(), LEFT), block_height/2, left_edge, thisIndex, rand());
          else
            thisProxy(receiver).receiveRefGhosts(iterations, map_nbr(thisIndex.getQuadI(), LEFT), block_height/2, left_edge);

          VB(logFile << ", " << map_nbr(thisIndex.getQuadI(), LEFT) << std::endl;);
          }
      else if(dir==RIGHT){
        for(int j=1; j<=block_height; j+=2) {
          right_edge[j/2] = (u[index(block_width-1,j)] + u[index(block_width,j)] + u[index(block_width-1,j+1)] +u[index(block_width, j+1)])/4;
          VB(logFile << right_edge[j/2] << "\t";);
            }
        VB(logFile << std::endl;);
        VB(logFile <<  map_nbr(thisIndex.getQuadI(), RIGHT) << std::endl;);
          if(!which)
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getQuadI(), RIGHT), block_height/2, right_edge, thisIndex, rand());
          else
            thisProxy(receiver).receiveRefGhosts(iterations, map_nbr(thisIndex.getQuadI(), RIGHT), block_height/2, right_edge);

          VB(logFile << ", " << map_nbr(thisIndex.getQuadI(), RIGHT) << std::endl;);
          }
      else if(dir==UP){
        for(int i=1; i<=block_width; i+=2){
          top_edge[i/2] = (u[index(i,1)] + u[index(i,2)] + u[index(i+1,1)] + u[index(i+1,2)])/4;
          VB(logFile << top_edge[i/2] << "\t";);
            }
        VB(logFile << std::endl;);
          if(!which)
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getQuadI(), UP), block_width/2, top_edge, thisIndex, rand());
          else
            thisProxy(receiver).receiveRefGhosts(iterations, map_nbr(thisIndex.getQuadI(), UP), block_width/2, top_edge);
          VB(logFile << ", " << map_nbr(thisIndex.getQuadI(), UP) << std::endl;);
          }
      else if(dir==DOWN){
        for(int i=1; i<=block_width; i+=2){
          bottom_edge[i/2] = (u[index(i,block_height-1)] + u[index(i,block_height)] + u[index(i+1,block_height-1)] + u[index(i+1,block_height)])/4;
          VB(logFile << bottom_edge[i/2] << "\t";);
            }
        VB(logFile << std::endl;);
          if(!which)
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getQuadI(), DOWN), block_width/2, bottom_edge, thisIndex, rand());
          else
            thisProxy(receiver).receiveRefGhosts(iterations, map_nbr(thisIndex.getQuadI(), DOWN), block_width/2, bottom_edge);

          VB(logFile << ", " << map_nbr(thisIndex.getQuadI(), DOWN) << std::endl;);
          }
      else {
        CkPrintf("noop send\n");
        return false;
      }
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

  for(int i=1; i<=block_width; i++)
    for(int j=1; j<=block_height; j++)
      u2[index(i,j)] = u[index(i,j)];

  for(int i=1; i<=block_width; i++)
    for(int j=1; j<=block_height; j++)
      u3[index(i,j)] = u[index(i,j)];


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

void Advection::process(int iter, int dir, int size, double gh[]){
//printf("[%d] process %d %d\n", thisIndex, iter, dir);
  VB(logFile << thisIndex.getIndexString() << " received data for direction " << dir << ", iteration " << iter << ", " << iterations << std::endl;);
    for(int i=0; i<size; i++){
      VB(logFile << gh[i] << '\t';);
        }
    VB(logFile << std::endl;);
    switch(dir){
    case LEFT:
      imsg++;
      VB(logFile << "Received From Left" << std::endl;);
        hasReceived.insert(LEFT);
      for(int i=0; i<size; i++)
        u[index(0,i+1)] = gh[i];
      //u[i+1][0] = gh[i];
      break;

    case RIGHT:
      imsg++;
      VB(logFile << "Received From RIGHT" << std::endl;);
        hasReceived.insert(RIGHT);
      for(int i=0; i<size; i++)
        u[index(block_width+1,i+1)] = gh[i];
      //u[i+1][block_width+1]=gh[i];
      break;

    case UP:
      imsg++;
      VB(logFile << "Received From UP" << std::endl;);
        hasReceived.insert(UP);
      for(int i=0; i<size; i++){
        u[index(i+1,0)] = gh[i];
        //u[block_height+1][i+1]=gh[i];
        VB(logFile << gh[i] << "\t" << std::endl;);
          }
      break;

    case DOWN:
      imsg++;
      VB(logFile << "Received From Down" << std::endl;);
        hasReceived.insert(DOWN);
      for(int i=0; i<size; i++)
        u[index(i+1,block_height+1)] = gh[i];
      //u[0][i+1]=gh[i];
      break;

    case LEFT_UP:
      imsg+=0.5;
      hasReceived.insert(LEFT_UP);
      for(int i=0; i<size; i++)
        u[index(0,i+1)] = gh[i];
            
      break;

    case LEFT_DOWN:
      imsg+=0.5;
      hasReceived.insert(LEFT_DOWN);
      for(int i=0; i<size; i++)
        u[index(0,block_height/2+i+1)] = gh[i];
      break;

    case RIGHT_UP:
      imsg+=0.5;
      hasReceived.insert(RIGHT_UP);
      for(int i=0; i<size; i++)
        u[index(block_width+1, i+1)]=gh[i];
      break;

    case RIGHT_DOWN:
      imsg+=0.5;
      hasReceived.insert(RIGHT_DOWN);
      for(int i=0; i<size; i++)
        u[index(block_width+1, block_height/2+i+1)]=gh[i];
      break;

    case UP_LEFT:
      imsg+=0.5;
      hasReceived.insert(UP_LEFT);
      for(int i=0; i<size; i++)
        u[index(i+1, 0)]=gh[i];
      break;

    case UP_RIGHT:
      imsg+=0.5;
      hasReceived.insert(UP_RIGHT);
      for(int i=0; i<size; i++)
        u[index(block_width/2+i+1,0)]=gh[i];
      break;

    case DOWN_LEFT:
      imsg+=0.5;
      hasReceived.insert(DOWN_LEFT);
      for(int i=0; i<size; i++)
        u[index(i+1,block_height+1)]=gh[i];
      break;

    case DOWN_RIGHT:
      imsg+=0.5;
      hasReceived.insert(DOWN_RIGHT);
      for(int i=0; i<size; i++)
        u[index(block_width/2+i+1,block_height+1)]=gh[i];
      break;

    default:
      CkAbort("ERROR\n");
    }
}

void Advection::sendReadyData2RefiningNeighbors(){
  for(int i=NUM_NEIGHBORS; i<3*NUM_NEIGHBORS; i++) {
    bool shouldSend = false;
    if(nbr_decision[i]==REFINE && !nbr_dataSent[i]) {
      if(i==RIGHT_UP){
        if(hasReceived.find(RIGHT_UP)!=hasReceived.end() &&
           (hasReceived.find(UP_RIGHT)!=hasReceived.end() || hasReceived.find(UP)!=hasReceived.end())){
          shouldSend = true;
        }
      }
      else if(i==RIGHT_DOWN){
        if(hasReceived.find(RIGHT_DOWN)!=hasReceived.end() &&
           (hasReceived.find(DOWN_RIGHT)!=hasReceived.end() || hasReceived.find(DOWN)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if (i==DOWN_RIGHT){
        if(hasReceived.find(DOWN_RIGHT)!=hasReceived.end() &&
           (hasReceived.find(RIGHT_DOWN)!=hasReceived.end() || hasReceived.find(RIGHT)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if(i==DOWN_LEFT){
        if(hasReceived.find(DOWN_LEFT)!=hasReceived.end() &&
           (hasReceived.find(LEFT_DOWN)!=hasReceived.end() || hasReceived.find(LEFT)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if(i==LEFT_DOWN){
        if(hasReceived.find(LEFT_DOWN)!=hasReceived.end() &&
           (hasReceived.find(DOWN_LEFT)!=hasReceived.end() || hasReceived.find(DOWN)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if(i==LEFT_UP){
        if(hasReceived.find(LEFT_UP)!=hasReceived.end() &&
           (hasReceived.find(UP_LEFT)!=hasReceived.end() || hasReceived.find(UP)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if(i==UP_LEFT){
        if(hasReceived.find(UP_LEFT)!=hasReceived.end() &&
           (hasReceived.find(LEFT_UP)!=hasReceived.end() || hasReceived.find(LEFT)!=hasReceived.end())){
          shouldSend = true;
        }
      }else if(i==UP_RIGHT){
        if(hasReceived.find(UP_RIGHT)!=hasReceived.end() &&
           (hasReceived.find(RIGHT_UP)!=hasReceived.end() || hasReceived.find(RIGHT)!=hasReceived.end())){
          shouldSend = true;
        }
      }
      if (shouldSend) {
        interpolateAndSend(i);
        nbr_dataSent[i]=true;
        ////terminator->msgSent(getRefinedNeighbor(i));
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
      bool shouldSend = false;
      if(i==RIGHT){
        if( hasReceived.find(RIGHT_UP)!=hasReceived.end() &&
            hasReceived.find(RIGHT_DOWN)!=hasReceived.end() &&
            (hasReceived.find(UP_RIGHT)!=hasReceived.end() || hasReceived.find(UP)!=hasReceived.end()) &&
            (hasReceived.find(DOWN_RIGHT) != hasReceived.end() || hasReceived.find(DOWN)!= hasReceived.end())
          )
        {
          shouldSend = true;
        }
      }
      else if(i==LEFT){
        if(hasReceived.find(LEFT_UP)!=hasReceived.end() &&
           hasReceived.find(LEFT_DOWN)!=hasReceived.end() &&
           (hasReceived.find(UP_LEFT)!=hasReceived.end() || hasReceived.find(UP)!=hasReceived.end()) &&
           (hasReceived.find(DOWN_LEFT) != hasReceived.end() || hasReceived.find(DOWN)!= hasReceived.end())
          ){
          shouldSend = true;
        }
      }
      else if(i==UP){
        if(hasReceived.find(UP_LEFT)!=hasReceived.end() &&
           hasReceived.find(UP_RIGHT)!=hasReceived.end() &&
           (hasReceived.find(RIGHT_UP)!=hasReceived.end() || hasReceived.find(RIGHT)!=hasReceived.end()) &&
           (hasReceived.find(LEFT_UP)!=hasReceived.end() || hasReceived.find(LEFT)!=hasReceived.end())
          ){
          shouldSend = true;
        }
      }
      else{ // if i==DOWN
        if(hasReceived.find(DOWN_LEFT)!=hasReceived.end() &&
           hasReceived.find(DOWN_RIGHT)!=hasReceived.end() &&
           (hasReceived.find(RIGHT_DOWN)!=hasReceived.end() || hasReceived.find(RIGHT)!=hasReceived.end()) &&
           (hasReceived.find(LEFT_DOWN)!=hasReceived.end() || hasReceived.find(LEFT)!=hasReceived.end())
          ){
          shouldSend = true;
        }
      }
      if (shouldSend) {
        interpolateAndSend(i);
        nbr_dataSent[i]=true;
      }
    }
  /*if(imsg==4)
    compute_and_iterate();*/
}

QuadIndex Advection::getRefinedNeighbor(int NBR) {
  QuadIndex receiver;
  if(NBR==DOWN_LEFT){
    receiver = nbr[DOWN].getChild(map_child(UP_LEFT));
  } else if(NBR==DOWN_RIGHT) {
    receiver = nbr[DOWN].getChild(map_child(UP_RIGHT));
  } else if(NBR==RIGHT_UP) {
    receiver = nbr[RIGHT].getChild(map_child(LEFT_UP));
  } else if(NBR==RIGHT_DOWN){
    receiver = nbr[RIGHT].getChild(map_child(LEFT_DOWN));
  } else if(NBR==UP_LEFT){
    receiver = nbr[UP].getChild(map_child(DOWN_LEFT));
  } else if(NBR==UP_RIGHT){
    receiver = nbr[UP].getChild(map_child(DOWN_RIGHT));
  } else if(NBR==LEFT_UP){
    receiver = nbr[LEFT].getChild(map_child(RIGHT_UP));
  } else if(NBR==LEFT_DOWN){
    receiver = nbr[LEFT].getChild(map_child(RIGHT_DOWN));
  } else {
    CkAbort("This only applies to refined neighbors");
  }
  return receiver;
}
    
void Advection::interpolateAndSend(int NBR){
  double sx, sy;
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

  QuadIndex receiver = getRefinedNeighbor(NBR);

  if(NBR==DOWN_LEFT){
    for(int i=0; i<block_width/2; i++){
      sx = (u[index(i,block_height)]-u[index(i+2, block_height)])/(2*dx);
      sy = (u[index(i+1, block_height-1)] - u[index(i+1, block_height+1)])/(2*dy);

      bottom_edge[wrap(2*i, block_width)] = u[index(i+1,block_height)] + sx*(dx/4) - sy*(dy/4);
      bottom_edge[wrap(2*i+1, block_width)] = u[index(i+1, block_height)] - sx*(dx/4) - sy*(dy/4);
    }
    receiver = nbr[DOWN].getChild(map_child(UP_LEFT));
    int val=rand();
    thisProxy(receiver).receiveGhosts(iterations, UP, block_width, bottom_edge,
                                      thisIndex, rand());
    VB(logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;);
#ifdef LOGGER
      for(int i=0; i<block_width; i++){
        logFile << bottom_edge[i] << '\t';
      }
    logFile << std::endl;
#endif
  }
  else if(NBR==DOWN_RIGHT){
    for(int i=block_width/2; i< block_width; i++){
      sx = (u[index(i,block_height)]-u[index(i+2, block_height)])/(2*dx);
      sy = (u[index(i+1, block_height-1)] - u[index(i+1, block_height+1)])/(2*dy);

      bottom_edge[wrap(2*i, block_width)] = u[index(i+1,block_height)] + sx*(dx/4) - sy*(dy/4);
      bottom_edge[wrap(2*i+1, block_width)] = u[index(i+1, block_height)] - sx*(dx/4) - sy*(dy/4);
    }

    receiver = nbr[DOWN].getChild(map_child(UP_RIGHT));
    thisProxy(receiver).receiveGhosts(iterations, UP, block_width, bottom_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_width; i++){
      logFile << bottom_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }
  else if(NBR==RIGHT_UP){
#ifdef LOGGER
    for(int j=0; j<block_height+2; j++){
      for(int i=0; i<3; i++)
        logFile << u[index(block_width-1+i, j)] << "\t";
      logFile << std::endl;
    }
#endif
    for(int i=0; i<block_height/2; i++){
      sx = (u[index(block_width-1, i+1)] - u[index(block_width+1, i+1)])/(2*dx);
      sy = (u[index(block_width, i)] - u[index(block_width, i+2)])/(2*dy);


      right_edge[wrap(2*i, block_height)] = u[index(block_width, i+1)] - sx*(dx/4) + sy*(dy/4);
      right_edge[wrap(2*i+1, block_height)] = u[index(block_width, i+1)] - sx*(dx/4) - sy*(dy/4);
    }
    receiver = nbr[RIGHT].getChild(map_child(LEFT_UP));//LEFT_UP child of the nieghbor
    thisProxy(receiver).receiveGhosts(iterations, LEFT, block_height, right_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_height; i++){
      logFile << right_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }
  else if(NBR==RIGHT_DOWN){
    for(int i=block_height/2; i<block_height; i++){
      sx = (u[index(block_width-1, i+1)] - u[index(block_width+1, i+1)])/(2*dx);
      sy = (u[index(block_width, i)] - u[index(block_width, i+2)])/(2*dy);

      right_edge[wrap(2*i, block_height)] = u[index(block_width, i+1)] - sx*(dx/4) + sy*(dy/4);
      right_edge[wrap(2*i+1, block_height)] = u[index(block_width, i+1)] - sx*(dx/4) - sy*(dy/4);
    }
    receiver = nbr[RIGHT].getChild(map_child(LEFT_DOWN));//LEFT_DOWN child of the neighbor
    thisProxy(receiver).receiveGhosts(iterations, LEFT, block_height, right_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_height; i++){
      logFile << right_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }
  else if(NBR==UP_LEFT){
    for(int i=0; i<block_width/2; i++){
      sx = (u[index(i,1)]-u[index(i+2,1)])/(2*dx);
      sy = (u[index(i+1,0)]-u[index(i+1,2)])/(2*dy);

      top_edge[wrap(2*i, block_width)] = u[index(i+1, 1)] + sx*(dx/4) + sy*(dy/4);
      top_edge[wrap(2*i+1, block_width)] = u[index(i+1, 1)] - sx*(dx/4) + sy*(dy/4);
    }
    receiver = nbr[UP].getChild(map_child(DOWN_LEFT));
    thisProxy(receiver).receiveGhosts(iterations, DOWN, block_width, top_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_width; i++){
      logFile << top_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }
  else if(NBR==UP_RIGHT){
    for(int i=block_width/2; i<block_width;i++){
      sx = (u[index(i,1)]-u[index(i+2,1)])/(2*dx);
      sy = (u[index(i+1,0)]-u[index(i+1,2)])/(2*dy);

      top_edge[wrap(2*i, block_width)] = u[index(i+1, 1)] + sx*(dx/4) + sy*(dy/4);
      top_edge[wrap(2*i+1, block_width)] = u[index(i+1, 1)] - sx*(dx/4) + sy*(dy/4);
    }
    receiver = nbr[UP].getChild(map_child(DOWN_RIGHT));
    thisProxy(receiver).receiveGhosts(iterations, DOWN, block_width, top_edge, thisIndex, rand());
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_width; i++){
      logFile << top_edge[i] << '\t';
    }
    logFile << std::endl;
  }
  else if(NBR==LEFT_UP){
    for(int i=0; i<block_height/2; i++){
      sx = (u[index(0, i+1)] - u[index(2,i+1)])/(2*dx);
      sy = (u[index(1, i)] - u[index(1, i+2)])/(2*dy);

      left_edge[wrap(2*i, block_height)] = u[index(1,i+1)] + sx*(dx/4) + sy*(dy/4);
      left_edge[wrap(2*i+1, block_height)] = u[index(1,i+1)] + sx*(dx/4) - sy*(dy/4);
    }
    receiver = nbr[LEFT].getChild(map_child(RIGHT_UP));//LEFT_UP child of the nieghbor
    thisProxy(receiver).receiveGhosts(iterations, RIGHT, block_height, left_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_height; i++){
      logFile << left_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }
  else if(NBR==LEFT_DOWN){
    for(int i=block_height/2; i<block_height; i++){
      sx = (u[index(0, i+1)] - u[index(2,i+1)])/(2*dx);
      sy = (u[index(1, i)] - u[index(1, i+2)])/(2*dy);

      left_edge[wrap(2*i, block_height)] = u[index(1,i+1)] + sx*(dx/4) + sy*(dy/4);
      left_edge[wrap(2*i+1, block_height)] = u[index(1,i+1)] + sx*(dx/4) - sy*(dy/4);
    }
    receiver = nbr[LEFT].getChild(map_child(RIGHT_DOWN));//LEFT_DOWN child of the neighbor
    thisProxy(receiver).receiveGhosts(iterations, RIGHT, block_height, left_edge, thisIndex, rand());
#ifdef LOGGER
    logFile << thisIndex.getIndexString() << " sending interpolated data to " << receiver.getIndexString() << std::endl;
    for(int i=0; i<block_height; i++){
      logFile << left_edge[i] << '\t';
    }
    logFile << std::endl;
#endif
  }

  ////terminator->msgSent(receiver);
};

void Advection::compute_and_iterate(){
  //logFile << "dt: " << dt << " ap:" << ap << " an:" << an << std::endl;
  //logFile << iterations << std::endl;
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
      /*This computation phase can be tested for correctness by running 
        extreme cases - like everyone wants to refine, 
        everyone wants to derefine, nobody wants to do anything */
      VB(logFile << "Entering Mesh Restructure Phase on " << thisIndex.getIndexString() << ", iteration " << iterations << std::endl;);
//      CkPrintf("%s in long phase 0 iteration %d\n", thisIndex.getIndexString().c_str(), iterations);
      contribute(CkCallback(CkIndex_Advection::startRemesh(), thisProxy));
      //startRemesh();
    }
    else {
      doStep();
    }
  }
  else {
    CkPrintf("Contribute from %s\n", thisIndex.getIndexString().c_str());

    //contribute();
  }
}

DECISION Advection::getGranularityDecision(){
  /*if(strcmp(thisIndex.getIndexString(),"00")==0 || strcmp(thisIndex.getIndexString(),"01")==0)// && iterations <15)
    return REFINE;
    else return STAY;*/
  //return REFINE;
  /*if(strcmp(thisIndex.getIndexString(),"00")==0)
    return REFINE;
    if(strlen(thisIndex.getIndexString())==4)
    return REFINE;
    if(iterations%50==0){
    return DEREFINE;
    }
    return STAY;*/
  /*if(iterations==5 || iterations==10 || iterations==15)
    return REFINE;
    return STAY;*/
    
  //if(thisIndex.getDepth()==min_depth)
  //  return STAY;
    
  /*if(strcmp(thisIndex.getIndexString(),"0111")==0)
    return REFINE;
    else return STAY;*/
    
  if(xc <= 0.10){
    return REFINE;
    }
    else if (xc <= 0.66 || thisIndex.getDepth()==min_depth){
    return STAY;
    }
    else return DEREFINE;

  for(int j=1; j<=block_height; j++)
    for(int i=1; i<=block_width; i++){
      if(u[index(i,j)]-u[index(i-1,j)]>0.33)
        return REFINE;
    }
  return STAY;

  if(rand()%5==0){
    return REFINE;
  }
  else{
    if(rand()%2==0 && thisIndex.getDepth()!=min_depth)
      return DEREFINE;	
    else 
      return STAY;
  }
/*
  double delx, dely, dely_f, ndim=2, ndim2=ndim*ndim;
  delx = 0.5*dx;
  dely = 0.5*dy;
  double refine_filter = 0.01;
  dely_f = dely;
  double error=0;
  double refine_cutoff=0.8, derfine_cutoff=0.2;

  for(int i=1; i <= block_width; i++){
  for(int j=1; j<=block_height; j++){
  delu[1][i][j] = u[i+1][j] - u[i-1][j];
  delu[1][i][j] = delu[1][i][j]*delx;

  delua[1][i][j] = abs(u[i+1][j]) + abs(u[i-1][j]);
  delua[1][i][j] = delua[1][i][j]*delx;

  // d/dy
  delu[2][i][j] = u[i][j+1] - u[i][j-1];
  delu[2][i][j] = delu[2][i][j]*dely_f;

  delua[2][i][j] = abs(u[i][j+1]) + abs(u[i][j-1]);
  delua[2][i][j] = delua[2][i][j]*dely_f;

  }
  }
    
  int istart=2, iend=block_width-1, jstart=2, jend=block_height-1;
  for (i=istart;i<=iend;i++){
  for (j=jstart;j<=jend;j++){

  delu2[1] = delu[1][i+1][j] - delu[1][i-1][j];
  delu2[1] = delu2[1]*delx;

  delu3[1] = abs(delux[1][i+1][j]) + abs(delux[1][i-1][j]);
  delu3[1] = delu3[1]*delx;

  delu4[1] = delua[1][i+1][j] + delua[1][i-1][j];
  delu4[1] = delu4[1]*delx;

  // d/dydx
  delu2[2] = delu[1][i][j+1] - delu[1][i][j-1];
  delu2[2] = delu2[2]*dely_f;

  delu3[2] = abs(delu[1][i][j+1]) + abs(delu[1][i][j-1]);
  delu3[2] = delu3[2]*dely_f;

  delu4[2] = delua[1][i][j+1] + delua[1][i][j-1];
  delu4[2] = delu4[2]*dely_f;

  // d/dxdy
  delu2[3] = delu[2][i+1][j] - delu[2][i-1][j]
  delu2[3] = delu2[3]*delx;

  delu3[3] = abs(delu[2][i+1][j]) + abs(delu[2][i-1][j]);
  delu3[3] = delu3[3]*delx;

  delu4[3] = delua[2][i+1][j] + delua[2][i-1][j];
  delu4[3] = delu4[3]*delx;

  // d/dydy
  delu2[4] = delu[2][i][j+1] - delux[2][i][j-1];
  delu2[4] = delu2[4]*dely_f;

  delu3[4] = abs(delu[2][i][j+1]) + abs(delu[2][i][j-1]);
  delu3[4] = delu3[4]*dely_f;

  delu4[4] = delua[2][i][j+1] + delua[2][i][j-1];
  delu4[4] = delu4[4]*dely_f;
		
  // compute the error

  double num = 0.;
  double denom = 0.;

  for (int kk=1;kk<= ndim2; kk++){  // kk= 1, 2, 3, 4
  num = num + pow(delu2[kk],2.);
  denom = denom + pow(delu3[kk]+ (refine_filter*delu4[kk]),2);
  // refine_filter is the epsilon in my writing sheet
  // refine_filter = 0.01 by default
  }
  // compare the square of the error
  if (denom == 0. && num != 0.){
  error = std::numeric_limits<double>::max(); // not sure this fn exist in C
					      // which means the largest number that is not infinity
                                              }else if (denom != 0.0){
                                              error = std::max(error, num/denom)
                                              }

                                              }
                                              }

                                              error = sqrt(error)
                                              if(error < derefine_cutoff && thisIndex.getDepth() > min_depth){
                                              return DEREFINE;
                                              }
                                              else if(error > refine_cutoff && thisIndex.getDepth() < max_depth){
                                              return REFINE;  
                                              }
                                              else return STAY;
*/

}

void Advection::resetMeshRestructureData(){
  /*Reset old data*/
  /*Phase1 Resetting*/
  for(int i=0; i<3*NUM_NEIGHBORS; i++)
    nbr_decision[i]=INV;

  decision = INV;

  for(int i=0; i<NUM_CHILDREN; i++)
    child_decision[i]=INV;
  VB(logFile << "setting parentHasAlreadyMadeDecision to false" << std::endl;);
    parentHasAlreadyMadeDecision=false;
  hasReceivedParentDecision=false;
  hasCommunicatedSTAY=false;
  hasCommunicatedREFINE=false;

  /*Phase2 resetting*/
  hasInitiatedPhase2 = false;
  hasAllocatedMemory=false;
}

void Advection::doMeshRestructure(){
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
    if(decision==INV)
      decision = getGranularityDecision();
    else if(decision==REFINE);
    else{//decision == STAY
      DECISION dec = getGranularityDecision();
      if(dec==STAY || dec==REFINE)
        decision=dec;
    }
    //ckout << thisIndex.getIndexString().c_str() << " decision = " << decision << ", iteration = " << iterations << endl;
    //initiate Phase1 of the computation
    if(decision==REFINE && !hasCommunicatedREFINE){
      hasCommunicatedREFINE=true;
      communicatePhase1Msgs();
    }
    else if(decision==STAY && !hasCommunicatedSTAY){
      hasCommunicatedSTAY=true;
      communicatePhase1Msgs();
    }
    //Inform Parent About Start of the Restructure Phase 
    // Think of a better way to do this because right now as many as 
    // number of unrefined children are informing the parent about the
    // start of the restructure phase while only one of them need to do so
    if(parent!=thisIndex) {
      //CkPrintf("%s -> %s startRemesh %d\n", thisIndex.getIndexString().c_str(),
      //parent.getIndexString().c_str(), iterations);
      //thisProxy(parent).startRemesh();
      ////terminator->msgSent(parent);
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

void Advection::phase1DoneQD() {
  ckout << thisIndex.getIndexString().c_str() << " signaled by QD" << endl;
  //contribute(CkCallback(CkIndex_Advection::doPhase2(), thisProxy));
}

/***** PHASE1 FUNCTIONS****/
void Advection::communicatePhase1Msgs(){
    
    
  if(decision==REFINE||decision==STAY){
    //tell the neighbor if it exists, also tell to children of the neighbor 
    // if that neighbor is refined
    //In case the neighbor does not exist, just tell to the parent of the 
    //non-existing neighbor

    for(int i=0; i<NUM_NEIGHBORS; i++){
      if(nbr_exists[i] && !nbr_isRefined[i]){
        VB(logFile << thisIndex.getIndexString() << " sending decision " << decision << " to " << nbr[i].getIndexString() << std::endl;);
          thisProxy(nbr[i]).exchangePhase1Msg(SENDER_DIR[i],decision);//Since Phase1Msgs are only refinement messages
          //terminator->msgSent(nbr[i]);
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
        //terminator->msgSent(q1);
        //terminator->msgSent(q2);
      }
      else{//send to the parent of the non-existing neighbor
        VB(logFile << thisIndex.getIndexString() << " sending decision " << decision << " to " << nbr[i].getParent().getIndexString() << std::endl;);
        thisProxy(nbr[i].getParent()).exchangePhase1Msg(map_nbr(thisIndex.getQuadI(), i), decision);
        //terminator->msgSent(nbr[i].getParent());
      }
    }
  }
  else{//No Need to do anything, just wait for 
    //neighbors to tell if they wish to derefine
  }

  //If my DECISION is to stay or to REFINE, tell the parent
  if((decision == REFINE || decision == STAY) && (parent!=thisIndex)){
    thisProxy(parent).informParent(thisIndex.getChildNum(), decision);
    //terminator->msgSent(parent);
  }
}

void Advection::informParent(int childNum, DECISION dec){//Will be called from two contexts: 
  //a) If the parent is a grandparent and also have children that are leaves
  //b) when a child sends REFINE/STAY message to the parent
  if(!hasReset){
    hasReset=true;
    resetMeshRestructureData();
  }
  VB(logFile << thisIndex.getIndexString() << ": in informParent called by child " << childNum << " with decision = " << dec << std::endl;);
    if(dec==REFINE){
      child_isRefined[childNum]=true;
      isGrandParent = true;
    }
  if(parentHasAlreadyMadeDecision==false){
    VB(logFile << "settin parentHasAlreadyMadeDecision to true " << std::endl;);
      parentHasAlreadyMadeDecision=true;
    //tell rest of the children which are not refined
    for(int i=0 ;i<NUM_CHILDREN; i++){
      if(i!=childNum && !child_isRefined[i]) {
        thisProxy(thisIndex.getChild(i)).recvParentDecision();
        //terminator->msgSent(thisIndex.getChild(i));
      }
    }
    //now tell the neighbors that I am Not Going to Derefine
    /*for(int i=0; i<NUM_NEIGHBORS; i++){//all the neighbors must be existing as I am a Parent
      CkAssert(nbr_exists[i]);
      thisProxy(nbr[i]).recvNeighborDecision(SENDER_DIR[i]);
      }*/
  }

  // Was called as an entry method, not as a local method
  if (dec != INV) {
    //terminator->doneSending();
    //terminator->msgProcessed();
  }
}

void Advection::recvParentDecision(){
  VB(logFile << thisIndex.getIndexString() << " has received decision from parent " << std::endl;);
    if(!hasReset){
      hasReset=true;
      resetMeshRestructureData();
    }

  hasReceivedParentDecision=true;
  if(decision==DEREFINE){
    decision=STAY;
    if(!hasCommunicatedSTAY){
      hasCommunicatedSTAY=true;
      communicatePhase1Msgs();
    }
  }

  //terminator->doneSending();
  //terminator->msgProcessed();
}

/*void Advection::recvNeighborDecision(DIR dir){
  if(!hasReset){
  hasReset=true;
  resetMeshRestructureData();
  }

  //This Message can be sent only by refined neighbors
  CkAssert(nbr_exists[dir] && nbr_isRefined[dir]);
  nbr_decision[dir]=STAY;

  if(!isRefined){
  if(decision==REFINE);//the message will be communicated elsewhere
  else if(decision==DEREFINE){//now I cannot derefine and tell parent about it
  decision=STAY;
  if(parent!=thisIndex && !hasReceivedParentDecision){
  thisProxy(parent).informParent(thisIndex.getChildNum(), STAY);//communicate only if the 
  //parent's decision has not already been received
  hasReceivedParentDecision=true;//because I know the parent's decision now
  }
  }
  else{
  //the message will be communicated elsewhere
  }
  }
  else if(isRefined){//then the parent already knows that it cannot derefine
  // however, tell the children that the neighbor is not going to derefine
  //
  int c1, c2;
  thisIndex.getSiblingInDirection(dir, c1, c2);
  thisProxy(thisIndex.getChild(c1)).recvStatusUpdateFromParent(dir);
  thisProxy(thisIndex.getChild(c2)).recvStatusUpdateFromParent(dir);
  }
  }*/

/*void Advection::recvStatusUpdateFromParent(int dir){
  hasReceivedStatusFromParent[dir]=true;
  if(nbr_decision[dir]!=REFINE){
  nbr_decision[dir]=STAY;
  }
  }*/

void Advection::exchangePhase1Msg(int dir, DECISION dec){//Phase1 Msgs are either REFINE or STAY messages
  VB(logFile << thisIndex.getIndexString() << " received decision " << dec << " from direction " << dir << std::endl; );
  if(!hasReset){
    hasReset=true;
    resetMeshRestructureData();
  }
  if(nbr_decision[dir]!=REFINE){
    VB(logFile << "setting decision of neighbor in dir " << dir << " to " << dec << std::endl;);
    nbr_decision[dir]=dec;
  }

    
  if(decision==DEREFINE || decision==STAY || decision==INV){//if decision was refine it would already have been 
    //communicated and the neighbors decision cannot change my decision now
    //Now check if my Decision Changes Because of this Message
    if(dir==LEFT||dir==RIGHT||dir==UP|dir==DOWN && dec==REFINE){
      decision=STAY;
      if(!hasCommunicatedSTAY){
        hasCommunicatedSTAY=true;
        communicatePhase1Msgs();
      }
    }
    else if(dir==LEFT_UP||dir==LEFT_DOWN||dir==RIGHT_UP||dir==RIGHT_DOWN||
	    dir==UP_LEFT||dir==UP_RIGHT||dir==DOWN_LEFT||dir==DOWN_RIGHT){
      if(dec==REFINE){// I am going to refine
        decision =REFINE;
        hasCommunicatedREFINE=true;
        VB(logFile << thisIndex.getIndexString() << " has changed its decision to REFINE" << std::endl;);
          communicatePhase1Msgs();
      }else if(dec==STAY){//dec can be either REFINE or STAY
        decision=dec;
                
        VB(logFile << thisIndex.getIndexString() << " decision\'s now is STAY" << std::endl;);
          if(decision==STAY && !hasCommunicatedSTAY){
            hasCommunicatedSTAY=true;
            communicatePhase1Msgs();
          }
          else if(decision==REFINE && !hasCommunicatedREFINE){
            hasCommunicatedREFINE=true;
            communicatePhase1Msgs();
          }
      }
    }
  }

  //terminator->doneSending();
  //terminator->msgProcessed();
}

#define index_c(i,j) (int)((j)*(block_width/2) + i)
ChildDataMsg::ChildDataMsg(int cnum, double mt, double mdt, int iter, double* u, bool* nbr_exists, bool* nbr_isRefined, DECISION* nbr_decision){
  //logFile << child_nbr_isRefined << std::endl;
  for(int i=1; i<= block_width; i+=2){
    for(int j=1; j<=block_height; j+=2){
      int idx = index_c(i/2, j/2);
      child_u[idx] = u[index(i,j)];
      child_u[idx] += u[index(i+1,j)];
      child_u[idx] += u[index(i, j+1)];
      child_u[idx] += u[index(i+1,j+1)];
      child_u[idx] /= 4;
      //logFile << child_u[idx] << "\t";
    }
    //logFile << std::endl;
  }
  childNum = cnum;
  iterations=iter;
  myt=mt;
  mydt=mdt;
  CmiMemoryCheck();//logFile << child_nbr_exists << ", " << child_nbr_isRefined << ", " << child_nbr_decision << std::endl;
  memcpy(child_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
  CmiMemoryCheck();//logFile << NUM_NEIGHBORS << std::endl;
  memcpy(child_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
  CmiMemoryCheck();
  memcpy(child_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);

}

/**** PHASE2 FUNCTIONS ****/
void Advection::doPhase2(){
  //cout << thisIndex.getIndexString() << " starting phase 2 " << iterations << std::endl;
  hasReceived.clear();//clear it up to track the ghosts layered required for restructure
  for(int i=0; i<NUM_NEIGHBORS; i++)
    nbr_dataSent[i]=false;

  imsg=0;

  VB(logFile << thisIndex.getIndexString() << " Entering Phase2 " << std::endl;);
  CmiMemoryCheck();
  if(isRefined){
    decision=INV;
  }
  VB(logFile << thisIndex.getIndexString() << " decision = " << decision << std::endl;);
  //Update the decision of your neighbors
  if(decision==STAY || decision==REFINE || decision==DEREFINE)
    for(int i=0; i<NUM_NEIGHBORS; i++) {
      if(nbr_decision[i]==INV){
        //CmiMemoryCheck();
        VB(logFile <<  thisIndex.getIndexString() << " has Not received Any Message From Nbr " << i << std::endl;);
        if((nbr_exists[i] && !nbr_isRefined[i])||!nbr_exists[i])
          nbr_decision[i]= DEREFINE;
        else {//when neighbor exists and is refined
          int d1, d2;
          if(i==RIGHT){d1=RIGHT_UP; d2=RIGHT_DOWN;}
          else if(i==LEFT){d1=LEFT_UP; d2=LEFT_DOWN;}
          else if(i==UP){d1=UP_LEFT; d2=UP_RIGHT;}
          else if(i==DOWN){d1=DOWN_LEFT; d2=DOWN_RIGHT;}
          // logFile << d1 << ": " << nbr_decision[d1] << ", " << d2 << ": " << nbr_decision[d2]<< std::endl;
          if(nbr_decision[d1]==INV){
            nbr_decision[d1]=DEREFINE;
            nbr_decision[d2]=DEREFINE;
          }
        }
      }
    }

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
          ////terminator->msgSent(lastSent)
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
          ////terminator->msgSent(thisIndex);
        }
        if(nbr_decision[d2]==REFINE) {
          thisProxy[thisIndex].getAndSendGhost();
          ////terminator->msgSent(thisIndex);
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
            ////terminator->msgSent(lastSent)
            ;
        }
      }
  }

  //CkPrintf("%s middle of phase 2 %d\n", thisIndex.getIndexString().c_str(), iterations);

  if(isRefined && !isGrandParent && !parentHasAlreadyMadeDecision){//I am a parent(whose None of the Children Are Refined) and has to derefine
    //Get Data From the Children and extrapolate it
  }
  else if(decision==DEREFINE && !isRefined){//send data to the parent
    VB(logFile << thisIndex.getIndexString() << " Sending Values to Parent" << std::endl;;);
      size_t sz = ((block_height)*(block_width))/4;
    ChildDataMsg *msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS) ChildDataMsg(thisIndex.getChildNum(), myt, mydt, iterations, u, nbr_exists, nbr_isRefined, nbr_decision);

    thisProxy(parent).recvChildData(msg);
    ////terminator->msgSent(parent);
    //deallocate all your memory and destroy yourself
    VB(logFile << "Destroying " << thisIndex.getIndexString() << std::endl;);
    shouldDestroy = true;
      VB(logFile << "Done Destroying " << thisIndex.getIndexString() << std::endl;);
      //this->~Advection();
      }
  else if(decision==REFINE){
    VB(logFile << "Refine called on " << thisIndex.getIndexString() << std::endl;);
    thisProxy[thisIndex].getGhostsAndRefine();
    ////terminator->msgSent(thisIndex);
  }
    
  if((decision==STAY) || (isRefined && !isGrandParent && !parentHasAlreadyMadeDecision)){
    VB(logFile << "calling doStep after QD on myself:" << isRefined << ", " << isGrandParent << ", " << parentHasAlreadyMadeDecision<< endl;);
    //CkStartQD(CkIndex_Advection::doStep(), &thishandle);
    isActive = true;
  } else {
    isActive = false;
  }

  //Update the Status of Your Neighbors
  if(decision == STAY){
    VB(logFile << "Phase2: " << thisIndex.getIndexString() << " updating the Status of Neighbors" << std::endl;);
      for(int i=0; i<NUM_NEIGHBORS; i++){
        if(!nbr_exists[i]){
          if(nbr_decision[i]==DEREFINE){
            VB(logFile << "ERROR(" << thisIndex.getIndexString() << "): doPhase(): Uncle(" << i << ") Cannot DEREFINE while I want to STAY" << std::endl;);
              CkExit();
          }
          else if(nbr_decision[i]==REFINE){
            nbr_exists[i]=true;
            nbr_isRefined[i]=false;
          }
          else if(nbr_decision[i]==STAY){
            nbr_exists[i]=false;
          }
        }
        else if(nbr_exists[i] && !nbr_isRefined[i]){
          if(nbr_decision[i]==REFINE)
            nbr_isRefined[i]=true;
          else if(nbr_decision[i]==STAY){
            nbr_exists[i]=true;
            nbr_isRefined[i]=false;
          }
          else{// the neighbor is going to get derefined
            nbr_exists[i]=false;
          }
        }
        else if(nbr_exists[i] && nbr_isRefined[i]){
          int d1, d2;
          if(i==RIGHT){d1=RIGHT_UP; d2=RIGHT_DOWN;}
          else if(i==LEFT){d1=LEFT_UP; d2=LEFT_DOWN;}
          else if(i==UP){d1=UP_LEFT; d2=UP_RIGHT;}
          else if(i==DOWN){d1=DOWN_LEFT; d2=DOWN_RIGHT;}
                
          if(nbr_decision[d1]==STAY){
            nbr_exists[i]=true; nbr_isRefined[i]=true;
          }else if(nbr_decision[d1]==DEREFINE){
            nbr_exists[i]=true; nbr_isRefined[i]=false;
          }
        }
      }
  }
  else if(decision==REFINE){//I will Now become Inactive and therefore I need not Store Neighbor Status
    isRefined=true;
    isGrandParent=false;
  }else if(isRefined && !isGrandParent && !parentHasAlreadyMadeDecision){// parent going to destroy its children
    isRefined = false;
  }
  //logFile << thisIndex.getIndexString() << " decision = " << decision << std::endl;
  VB(logFile << "setting parentHasAlreadyMadeDecision to false" << endl;);
  parentHasAlreadyMadeDecision=false;
  hasReset=false;

  if (decision != REFINE) {
    //CkPrintf("%s nonrefine phase 2b iteration %d\n", thisIndex.getIndexString().c_str(), iterations);
    //    fflush(stdout);
//     contribute(CkCallback(CkReductionTarget(Advection, phase2Done), thisProxy));
  } else {
    //CkPrintf("%s refine phase 2e %d\n", thisIndex.getIndexString().c_str(), iterations);
  }  
}

void Advection::recvChildData(ChildDataMsg *msg){
  VB(logFile << "Mem Check at Beginning of recvChildData" << std::endl;);
  VB(logFile << thisIndex.getIndexString() << " received data from Child " << msg->childNum << " for coarsening" << std::endl;);
    myt = msg->myt;
  mydt = msg->mydt;
  iterations = msg->iterations;

  if(!hasAllocatedMemory){
    hasAllocatedMemory=true;
    mem_allocate_all();
  }
  int st_i, end_i, st_j, end_j;
  if(msg->childNum==0){
    st_i=block_width/2+1; end_i=block_width;
    st_j=1; end_j=block_height/2;
  }
  else if(msg->childNum==1){
    st_i=1; end_i=block_width/2;
    st_j=1; end_j=block_height/2;
  }
  else if(msg->childNum==2){
    st_i=1; end_i=block_width/2;
    st_j=block_height/2+1; end_j=block_height;
  }
  else if(msg->childNum==3){
    st_i=block_width/2+1; end_i=block_width;
    st_j=block_height/2+1;end_j=block_height;
  }
  else{
    VB(logFile << "Error: recvChildData(ChildDataMsg*) received " << msg->childNum << "as ChildNum" <<std::endl;);
      CkExit();
  }
  VB(logFile << "Check Memory 1" << thisIndex.getIndexString() << std::endl;);

    int ctr=0;
  for(int j=st_j; j<=end_j; j++){
    for(int i=st_i; i<=end_i; i++){
      u[index(i,j)]=msg->child_u[ctr];
      VB(logFile << msg->child_u[ctr] << ", " << u[index(i,j)] << "\t";);
        ctr++;
    }
    VB(logFile << std::endl;);
      }

  //Update the Status of Your Neighbors based on Data Sent from the Children
  int c1, c2;
  if(msg->childNum==0){
    c1=RIGHT; c2=UP;
  }
  else if(msg->childNum==1){
    c1=LEFT; c2=UP;
  }
  else if(msg->childNum==2){
    c1=LEFT; c2=DOWN;
  }
  else if(msg->childNum==3){
    c1=RIGHT; c2=DOWN;
  }
    
  this->iterations=msg->iterations;
  this->myt=msg->myt;
  this->mydt=msg->mydt;
  setNbrStatus(c1, msg);
  setNbrStatus(c2, msg);
    
  ////terminator->msgProcessed();

  CmiMemoryCheck();
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
    liveVizDeposit(m,0,0,0,0,NULL, this);
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
  double sx, sy;
  int m=1, n=1;
  for(int i=xstart; i<=xend; i++){
    for(int j=ystart; j<=yend; j++){
      sx = (-u[index(i-1,j)]+u[index(i+1,j)])/(2*dx);
      sy = (-u[index(i,j-1)]+u[index(i,j+1)])/(2*dy);

      refined_u[index_l(2*(m-1), 2*(n-1))] = u[index(i,j)] - sx*(dx/4) - sy*(dy/4);
      refined_u[index_l(2*(m-1)+1, 2*(n-1))] = u[index(i,j)] + sx*(dx/4) - sy*(dy/4);
      refined_u[index_l(2*(m-1), 2*(n-1)+1)] = u[index(i,j)] - sx*(dx/4) + sy*(dy/4);
      refined_u[index_l(2*(m-1)+1, 2*(n-1)+1)] = u[index(i,j)] + sx*(dx/4) + sy*(dy/4);
      n++;
    }
    m++;
    n=1;
  }
}

void Advection::refine(){
  //Spawn The four children and give them the data
  //Assuming we already have the new boundary data
    
  //Interpolate the data and give it to the children when they are initialized
  // boundaries of the children will have to be sent by the neighbor
  VB(logFile << thisIndex.getIndexString() << " is refining" << std::endl;);
    double sx, sy, s;
  size_t sz = (block_width)*(block_height); 
  vector<double> refined_u(sz);

  QuadIndex child01 = thisIndex.getChild("01"), child00 = thisIndex.getChild("00"),
    child10 = thisIndex.getChild("10"), child11 = thisIndex.getChild("11");


  interpolate(u, refined_u, 1, block_width/2, 1, block_height/2);
  //initialize the child
  InitRefineMsg * msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, xmin, ymin+(ny*dy)/2, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
  thisProxy(child01).insert(msg);

  interpolate(u, refined_u, block_width/2+1, block_width, 1, block_height/2);
  //initialize the child
  msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, xmin+(nx*dx)/2, ymin+(ny*dy)/2, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
  thisProxy(child00).insert(msg);

  interpolate(u, refined_u, 1, block_width/2, block_height/2+1, block_height);
  //init the child
  msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, xmin, ymin, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
  thisProxy(child10).insert(msg);

  interpolate(u, refined_u, block_width/2+1, block_width, block_height/2+1, block_height);
  //init the child
  msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, xmin+(nx*dx)/2, ymin, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
  thisProxy(child11).insert(msg);
  thisProxy.doneInserting();

  // //terminator->msgSent(child00);
  // //terminator->msgSent(child01);
  // //terminator->msgSent(child10);
  // //terminator->msgSent(child11);
  //CkPrintf("%s parent phase 2b iteration %d\n", thisIndex.getIndexString().c_str(), iterations);
  //fflush(stdout);
//   contribute(CkCallback(CkReductionTarget(Advection, phase2Done), thisProxy));

  VB(logFile << thisIndex.getIndexString() << " done with refinement" << std::endl;);;
}

Advection::Advection(InitRefineMsg* msg)
  : AdvTerm(thisProxy, thisIndex, true), CBase_Advection()
{
  //ckout << thisIndex.getIndexString().c_str() << " created 2" << endl;
  __sdag_init();
  isActive = true;
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
      if(nbr[dir].getParent() == thisIndex.getParent()){//if parents are same
        nbr_exists[dir]=true;
        nbr_isRefined[dir]=false;
      }
      else{//when the parents are not the same
        //if I want to refine then no-one can stop me from doing so
        //so it is possible that my parents neighbor do not exist
        //at this moment but a notification has been sent that they
        //should be generated
        if(msg->parent_nbr_exists[dir]){
          if(msg->parent_nbr_isRefined[dir]){
            nbr_exists[dir]=true;
            //check for the decision of the neighbors sent by parent
            nbr_isRefined[dir]=false;
            VB(logFile << thisIndex.getQuadI() << std::endl;);
                    
            VB(logFile << thisIndex.getIndexString() << " nbr in direction " << getNbrDir(thisIndex.getQuadI(), dir) << " has decision = " << msg->parent_nbr_decision[getNbrDir(thisIndex.getQuadI(), dir)] << std::endl;);
              if(msg->parent_nbr_decision[getNbrDir(thisIndex.getQuadI(), dir)]==REFINE){
                nbr_isRefined[dir]=true;
                VB(logFile << thisIndex.getIndexString() << " setting " << dir << " to refined" << std::endl;);
                  }
          }
          else{
            nbr_exists[dir]=false;
            nbr_isRefined[dir]=false;
            if(msg->parent_nbr_decision[dir]==REFINE)
              nbr_exists[dir]=true;
          }
        }
        else{
          nbr_exists[dir]=false;
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
  iterations = msg->iterations;

  mem_allocate_all();

  for (int i = 0; i < NUM_CHILDREN; ++i)
    child_isRefined[i] = false;
    
  //x[i] and y[i] need not be initialized they are needed only for setting initial conditions
  //delete [] x;
  //delete [] y;
 	
  //Initialize u - For boundaries I have to wait for the neighbors
  //to send the values, rest of it can be initialized by the values 
  //received from the parent
  int ctr=0;
  for(int j=1; j<=block_height; j++)
    for(int i=1; i<=block_width; i++)
      u[index(i,j)]=msg->refined_u[ctr++];

#ifdef LOGGER
  logFile << "New Child Values" << std::endl;
  for(int i=0; i<block_height; i++){
    for(int j=0; j<block_width; j++)
      logFile << u[index(j+1,i+1)] << "\t";
    logFile << std::endl;
  }
#endif
  //delete the message
  delete msg;

  ////terminator->msgProcessed();

  //CkPrintf("%s child phase 2b iteration %d\n", thisIndex.getIndexString().c_str(), iterations);
  //fflush(stdout);

  //call Quiesence detection to begin next iteration*/
  //CkStartQD(CkIndex_Advection::doStep(), &thishandle);
}
#include "Advection.def.h"
