#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>
#include <boost/assign/list_of.hpp>
#include "boost/filesystem.hpp"

using namespace std;
using namespace boost::assign;

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

int min_depth, max_depth;

#define wrap_x(a) (((a)+num_chare_rows)%num_chare_rows)
#define wrap_y(a) (((a)+num_chare_cols)%num_chare_cols)


CkArrayID a;

int nframe;
double xmin, xmax, ymin, ymax;
double xctr, yctr, radius;
double dx, dy, v;
double ap, an;
double tmax, t, dt, cfl;
map<string, DIR> nbrDirectionMap;
map<DIR, DIR> reverse_dir_map;

char* decimal_to_binary_string(int num, int len){
    char* _ret = new char[len+1];
    int i=0;
    _ret[len]=0;
    while(num){
        _ret[len-i-1] = (num&1==1)?'1':'0';
        num = num>>1;
        i++;
    };
    for(; i<len; i++)
        _ret[len-i-1]='0';
    return _ret;
}

Main:: Main(CkArgMsg* m){

    initUtils();
    mainProxy = thisProxy;
    boost::filesystem::remove_all("out"); boost::filesystem::remove_all("log");
    boost::filesystem::create_directory("out"); boost::filesystem::create_directory("log");

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
     dx = (xmax - xmin)/double(array_width);//ckout << "dx: " << dx << endl;
     dy = (ymax - ymin)/double(array_height);//ckout << "dy: " << dy << endl;
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
                                                        CkNumPes(), array_width, array_height);
     CkCallback *cb = new CkCallback(CkCallback::ckExit);

     qtree = CProxy_Advection::ckNew();

     qtree.ckSetReductionClient(cb);
     //save the total number of worker chares we have in this simulation
     num_chares = num_chare_rows*num_chare_cols;
     int depth = (int)(log(num_chares)/log(4));
     min_depth = depth;
     max_depth = 100;

     QuadIndex qindex;
     for(int i=0; i < num_chares; i++){
        char* str = decimal_to_binary_string(i, 2*depth);
        ckout << str << endl;
        qindex = *new QuadIndex(str);
        qtree[qindex].insert(xmin, xmax, ymin, ymax);
     }
     qtree[qindex].doneInserting();
     
     for(int i=0; i < num_chares; i++){
        char* str = decimal_to_binary_string(i, 2*depth);
        qindex = *new QuadIndex(str);
        qtree[qindex].doStep();
     }
     
     /*queue<QuadIndex> q;
     q.push("");
       
     for(int i=0; i<depth; i++){
        int size = q.size();
        for(int j=0; j<size; j++){
            qindex = q.front(); q.pop();
            qtree(qindex).refine();
            qtree.doneInserting();
            for(int dir=0; dir<4; dir++)
                q.push(qindex.getChild(dir));
        }
    }
    CkVec<QuadIndex> v = *new CkVec<QuadIndex>();
    qtree[qindex].refine();
    v.push_back(qindex);
    v.push_back(qindex.getChild("00")); v.push_back(qindex.getChild("01"));
    v.push_back(qindex.getChild("10")); v.push_back(qindex.getChild("11"));

    qindex = qindex.getChild("00");
    qtree[qindex].refine();
    v.push_back(qindex.getChild("00")); v.push_back(qindex.getChild("01"));
    v.push_back(qindex.getChild("10")); v.push_back(qindex.getChild("11"));
    
    thisProxy.printTreeInformation(v);*/
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

void Main::initUtils(){

    //first two characters in the key are the source's index
    //and the last two characters are the target's index
    nbrDirectionMap.insert(std::pair<string, DIR>("0001", LEFT));
    nbrDirectionMap.insert(std::pair<string, DIR>("0011", DOWN));
    nbrDirectionMap.insert(std::pair<string, DIR>("0100", RIGHT));
    nbrDirectionMap.insert(std::pair<string, DIR>("0110", DOWN));
    nbrDirectionMap.insert(std::pair<string, DIR>("1001", UP));
    nbrDirectionMap.insert(std::pair<string, DIR>("1011", RIGHT));
    nbrDirectionMap.insert(std::pair<string, DIR>("1100", UP));
    nbrDirectionMap.insert(std::pair<string, DIR>("1110", LEFT));

    reverse_dir_map = map_list_of (RIGHT, LEFT) (LEFT, RIGHT) (UP, DOWN) (DOWN, UP);
}

void Main::printTreeInformation(CkVec<QuadIndex> list){
    for(int i=0; i<list.size(); i++){
        QuadIndex qindex = list[i];
        qtree[qindex].printState();
    }   
}
#include "Main.def.h"
