#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>

using namespace std;

#include "charm++.h"
#include "liveViz.h"
#include "Main.h"

#include <algorithm>

CProxy_Main mainProxy;

int array_height;
int array_width;

int num_chare_rows;
int num_chare_cols;

int block_height;
int block_width;

#define wrap_x(a) (((a)+num_chare_rows)%num_chare_rows)
#define wrap_y(a) (((a)+num_chare_cols)%num_chare_cols)


CkArrayID a;

int nframe;
double xmin, xmax, ymin, ymax;
double xctr, yctr, radius;
int nx, ny;
double dx, dy, v;
double ap, an;
double tmax, t, dt, cfl;

Main:: Main(CkArgMsg* m){
    mainProxy = thisProxy;
    iterations = 0;

    if(m->argc < 3){
        ckout << "Usage: " << m->argv[0] << "[array_size] [block_size]" << endl; 
        CkExit();
    }
           
    array_height = array_width =atoi(m->argv[1]);
    block_height = block_width = atoi(m->argv[2]);

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
     io_outnum = 0;
     io_tnext = 0;
     nframe = 20;
     io_tout = tmax/nframe;

     //call colormap
     nx = array_height;
     ny = array_width;
     dx = (xmax - xmin)/double(nx);//ckout << "dx: " << dx << endl;
     dy = (ymax - ymin)/double(ny);//ckout << "dy: " << dy << endl;
     //ckout << min(dx, dy) << endl;	
     dt = min(dx,dy)/v * cfl;
     if ((t + dt) >= tmax )
           dt = tmax - t;
     t = t+dt;
     //ckout << "dt: " << dt << endl;
     xctr = 0.3;
     yctr = 0.5;
     radius = 0.2;

     ap = std::max(v, 0.0);
     an = std::min(v, 0.0);
     //ckout << "ap: " <<ap << endl;
     //ckout << "an: " <<an << endl;

     /*****End Initialization **********/
     CkPrintf("Running Advection on %d processors with (%d,%d) elements\n",
                                                        CkNumPes(), nx, ny);
     CkCallback *cb = new CkCallback(CkCallback::ckExit);

     qtree = CProxy_Advection::ckNew();

     qtree.ckSetReductionClient(cb);
     //save the total number of worker chares we have in this simulation
     num_chares = num_chare_rows*num_chare_cols;
     int depth = (int)(log(num_chares)/log(4));
      
     QuadIndex qindex = *new QuadIndex("");
     qtree[qindex].insert(false, false, false);
    
     queue<QuadIndex> q;
     q.push("");
       
     /*for(int i=0; i<depth; i++){
        int size = q.size();
        for(int j=0; j<size; j++){
            qindex = q.front(); q.pop();
            qtree(qindex).refine();
            qtree.doneInserting();
            for(int dir=0; dir<4; dir++)
                q.push(qindex.getChild(dir));
        }
    }*/
    CkVec<QuadIndex> v = *new CkVec<QuadIndex>();
    qtree[qindex].refine();
    v.push_back(qindex);
    v.push_back(qindex.getChild("00")); v.push_back(qindex.getChild("01"));
    v.push_back(qindex.getChild("10")); v.push_back(qindex.getChild("11"));

    qindex = qindex.getChild("00");
    qtree[qindex].refine();
    v.push_back(qindex.getChild("00")); v.push_back(qindex.getChild("01"));
    v.push_back(qindex.getChild("10")); v.push_back(qindex.getChild("11"));
    
    thisProxy.printTreeInformation(v);
    //setup - liveViz
    CkCallback c(CkIndex_Advection::requestNextFrame(0), qtree);
    liveVizConfig cfg(liveVizConfig::pix_color, true);
    liveVizInit(cfg, a, c);
      
    /*int size = q.size();
    ckout << "Size of Queue " << size << endl;
    for(int i=0; i<size; i++){
        qindex = q.front(); q.pop();
        ckout << "Calling doStep for " << qindex.getIndexString() << endl;
        qtree[qindex].doStep();
    }*/
    /*qtree[qindex].doStep();*/
}

void Main::printTreeInformation(CkVec<QuadIndex> list){
    for(int i=0; i<list.size(); i++){
        QuadIndex qindex = list[i];
        qtree[qindex].printState();
    }   
}
#include "Main.def.h"
