#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>

using namespace std;

#include "charm++.h"
#include "QuadIndex.cc"
#include "liveViz.h"
#include "advection.decl.h"
#include "Main.h"

//#include <algorithm>

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

double max(double a, double b){
    return (a>b)?a:b;
}

double min(double a, double b){
    return (a<b)?a:b;
}

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

     ap = max(v, 0.0);
     an = min(v, 0.0);
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
        
    //setup - liveViz
    CkCallback c(CkIndex_Advection::requestNextFrame(0), qtree);
    liveVizConfig cfg(liveVizConfig::pix_color, true);
    liveVizInit(cfg, a, c);
      
    int size = q.size();
    for(int i=0; i<size; i++){
        qindex = q.front(); q.pop();
        ckout << "Calling doStep for " << qindex.getIndexString() << endl;
        qtree[qindex].doStep();
    }
}

template<class T>
T square(T x){
  	return x*x;	   
}

class Advection: public CBase_Advection{
Advection_SDAG_CODE
public:
    //tree information
    bool isdummy;
    bool hasRealChildren;
    bool hasDummyChildren;

    QuadIndex nbr[4], parent;
    int xc, yc;

    //data
    int imsg;

    double** u;
    double** u2;
    double** u3;
    double *x;
    double *y;

    double *left_edge;
    double *right_edge;

    int iterations;
    
    double up;
    double un;

    void mem_allocate(double* &p, int size){
        p = new double[size];
    }

    void mem_allocate_all(){
        mem_allocate(u, (block_width+2)*(block_height+2));
        mem_allocate(u2, (block_width+2)*(block_height+2));
        mem_allocate(u3, (block_width+2)*(block_height+2));

        mem_allocate(x, block_width+2);
        mem_allocate(y, block_height+2);

        mem_allocate(left_edge, block_height);
        mem_allocate(right_edge, block_height);
    }

    Advection(bool isdummy, bool hasRealChildren, bool hasDummyChildren){
        CBase_Advection();
        this->isdummy = isdummy;
        this->hasRealChildren = hasRealChildren;
        this->hasDummyChildren = hasDummyChildren;
        
        for(int dir=UP; dir<=RIGHT; ++dir)
            nbr[dir] = thisIndex.getNeighbor(dir);
        if(thisIndex.nbits!=0)
        parent = thisIndex.getParent();

        ckout << "Constructor for " << thisIndex.getIndexString() << " called" << endl;
        advection();
    }
    
    
    void refine(){	
        // should always be called on a node not having Real Chilren
        ckout << "refine caled for " << thisIndex.getIndexString() << endl;
        if(hasDummyChildren){               
            hasDummyChildren = false;
            hasRealChildren = true;

            for(int i=0; i<4; i++){
                thisProxy(thisIndex.getChild(i)).setReal();ckout << "Setting " << thisIndex.getChild("00").getIndexString() << " as real" << endl;
            }
        }
        else if(!hasRealChildren){
            // Create new real chilren
            hasDummyChildren = false;
            hasRealChildren = true;

            for(int i=0; i<4; i++){
                thisProxy(thisIndex.getChild(i)).insert(false, false, false); ckout << "Creating real child: " << thisIndex.getChild("00").getIndexString() << endl; 
            }
            thisProxy.doneInserting();
        }
        
        for(int i=0; i<4; i++){
            if(thisIndex!=nbr[i])
                thisProxy(nbr[i]).inform_nbr_of_refinement();
        }
    }
    
    void inform_nbr_of_refinement(){
        ckout << "inform_nbr called for " << thisIndex.getIndexString() << endl; 
        if(isdummy){
            //ask the parent to refine itself
            thisProxy(parent).refine();
        }    
        if(!hasDummyChildren && !hasRealChildren){
            hasDummyChildren = true;
            hasRealChildren = false;

            //create dummy children  
            for(int i=0; i<4; i++){
                thisProxy(thisIndex.getChild(i)).insert(true, false, false);ckout << "Creating dummy child: " << thisIndex.getChild("00").getIndexString() << endl;
            }
            thisProxy.doneInserting();
        }
    }
    
    void setDummy(){
        isdummy = true;
    }

    void setReal(){
        isdummy = false;
    }

    Advection(){
        advection();
    }

    void advection(){
        __sdag_init();
        usesAtSync = CmiTrue;
        //ckout << "In Advection: " << endl;
        mem_allocate_all();
        iterations=0;

        thisIndex.getCoordinates(xc, yc);

        for(int i=0; i<block_width+2; i++){
            x[i] = xmin + double(xc*block_width+i)*dx - 0.5*dx;
           // ckout << x[i] << endl;
        }
        //ckout << endl;

        for(int i=0; i<block_height+2; i++){
            y[i] = ymin + double(yc*block_height+i)*dy - 0.5*dy;
            //ckout << y[i] << endl;
        }
        //ckout << endl;

        double rsq;
        //ckout << "In Adfvection2" << endl;
        for(int i=0; i<block_height+2; i++){
            for(int j=0; j<block_width+2; j++){
                rsq = square(x[i] - xctr) + square(y[j] - yctr);
                if(rsq <= square(radius))
                    u[i][j] = 2;
                else u[i][j] = 1;
            }
        }
        //ckout << "In Advection 3: " <<u[1][1]<< endl;
        /*for(int i=0; i<block_height; i++){
            for(int j=0; j<block_width; j++)
                ckout << u[i][j] << "\t";
            ckout << endl;
        }
        ckout << endl;*/
    }

    Advection(CkMigrateMessage* m) {__sdag_init();}

    //added for array migration - see how 2D arrays can be packed
    void pup(PUP::er &p){
        ckout << "In PUP" << endl;
        CBase_Advection::pup(p);
       __sdag_pup(p);

	//p | u;
	//p | u2;
	//p | u3;
        p | imsg;
    }
    
    void free_mem2D(double** &x){
        for(int i=0; i<block_height+2; i++)
            delete [] x[i];
        delete [] x;
    }

    void free_mem1D(double* &x){
        delete [] x;
    }

    ~Advection(){
        ckout << "In Destructor" << endl;
        free_mem2D(u);
        free_mem2D(u2);
        free_mem2D(u3);

        free_mem1D(x);
        free_mem1D(y);
        free_mem1D(left_edge);
        free_mem1D(right_edge);
    }

    void begin_iteration(void) {
        //ckout << "String: " << thisIndex.getIndexString() << endl;

        //ckout << "u[1][1]: " << u[1][1] << endl;
        for(int i=1; i<=block_width; i++)
            for(int j=1; j<=block_height; j++){
                u2[i][j] = u[i][j];
            }

        for(int i=1; i<=block_width; i++)
            for(int j=1; j<=block_height; j++)
                u3[i][j]=u[i][j];

        iterations++;

        for(int i=0; i<block_height; i++){
            left_edge[i] = u[1][i+1];
            right_edge[i] = u[block_width][i+1];
        }
        //ckout << "Left Neighbor of " << thisIndex.getIndexString() << " is " << nbr[LEFT].getIndexString() << endl;
        //send my left edge
        thisProxy(nbr[LEFT]).receiveGhosts(iterations, RIGHT, block_height, left_edge);

        //ckout << "Right Neighbor of " << thisIndex.getIndexString() << " is " << nbr[RIGHT].getIndexString() << endl;
        //send my right edge
        thisProxy(nbr[RIGHT]).receiveGhosts(iterations, LEFT, block_height, right_edge);


        //ckout << "Top Neighbor of " << thisIndex.getIndexString() << " is " << nbr[UP].getIndexString() << endl;
        //send my top edge
        thisProxy(nbr[UP]).receiveGhosts(iterations, DOWN, block_width, &u[1][1]);

        //ckout << "Bottom Neighbor of " << thisIndex.getIndexString() << " is " << nbr[DOWN].getIndexString() << endl;
        //send my bottom edge
        thisProxy(nbr[DOWN]).receiveGhosts(iterations, UP, block_width, &u[block_height][1]);

    }
    
    void process(int dir, int size, double gh[]){
        switch(dir){
            case LEFT:
                for(int i=0; i<size; i++)
                    u[0][i+1] = gh[i];
            break;

            case RIGHT:
                for(int i=0; i<size; i++)
                    u[block_width+1][i+1] = gh[i];
            break;

            case UP:
                for(int i=0; i<size; i++)
                    u[i+1][block_height+1] = gh[i];
                break;

             case DOWN:
                for(int i=0; i<size; i++)
                    u[i+1][0] = gh[i];
                break;
            default:
                CkAbort("ERROR\n");
        }
    }
    
    void compute_and_iterate(){
        ckout << iterations << endl;
        ckout << "In compute and iterate" << endl;
        for(int i=1; i<=block_height; i++){
            for(int j=1; j<=block_width; j++){
                up = (u[i+1][j] - u[i][j])/dx;
                un = (u[i][j]-u[i-1][j])/dx;

                u2[i][j] = u[i][j] - dt* (ap*un + an*up);
                
        //        ckout << u2[i][j] << endl;
            }
        }

        for(int i=1; i<=block_height; i++)
            for(int j=1; j<=block_width ;j++){
                up = (u[i][j+1] - u[i][j])/dy;
                un = (u[i][j] - u[i][j-1])/dy;

               u3[i][j] = u[i][j] - dt*(ap*un + an*up);
            }

        for(int i=1; i<=block_height; i++)
            for(int j=1; j<=block_width; j++)
                u[i][j] = 0.5*(u2[i][j] + u3[i][j]);
        
        for(int i=1; i<=block_height; i++){
            for(int j=1; j<=block_width; j++)
                ckout << u[i][j] << "\t";
                ckout << endl;
        }
        ckout << endl;
        iterate();
    }
    void iterate() {
             t = t+dt;
             if(t < tmax){
                 dt = min(dx,dy)/v * cfl;
                 if ((t + dt) >= tmax )
                     dt = tmax - t;
                 doStep();    
                 ckout << "Doing another iteration" << endl;
             }
             else
                 contribute();
    }
    
    void requestNextFrame(liveVizRequestMsg *m){
        ckout << "In request New Frame: ********************************" << endl;

        int sy = xc*block_width;
        int sx = yc*block_height;
        int w= block_width, h = block_height;

        unsigned char *intensity = new unsigned char[3*w*h];
        if(isdummy || hasRealChildren){
            liveVizDeposit(m, 0, 0, 0, 0, intensity, this);
            delete [] intensity;
            return;
        }

        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
               // ckout << u[i+1][j+1] << endl;
                intensity[3*(i*w+j)+ 0] = 255;// red component
                intensity[3*(i*w+j)+ 1] = 255 - (u[i+1][j+1]-1)*255;// BLUE component
                intensity[3*(i*w+j)+ 2] = 255 - (u[i+1][j+1]-1)*255;// green component
            }
        }
        liveVizDeposit(m, sx,sy, w,h, intensity, this);
        delete[] intensity;
    }
};

#include "advection.def.h"
