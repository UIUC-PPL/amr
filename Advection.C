#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <queue>

using namespace std;

#include "charm++.h"
#include "liveViz.h"
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


void Advection::mem_allocate(double* &p, int size){
    p = new double[size];
}

void Advection::mem_allocate_all(){
    mem_allocate(u, (block_width+2)*(block_height+2));
    mem_allocate(u2, (block_width+2)*(block_height+2));
    mem_allocate(u3, (block_width+2)*(block_height+2));

    mem_allocate(x, block_width+2);
    mem_allocate(y, block_height+2);

    mem_allocate(left_edge, block_height);
    mem_allocate(right_edge, block_height);
}

Advection::Advection(bool isdummy, bool hasRealChildren, bool hasDummyChildren){
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
    
    
void Advection::refine(){	
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

void Advection::inform_nbr_of_refinement(){
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
    
void Advection::advection(){
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
            rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[i]-yctr);
            if(rsq <= radius*radius)
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

    //added for array migration - see how 2D arrays can be packed
void Advection::pup(PUP::er &p){
    ckout << "In PUP" << endl;
    CBase_Advection::pup(p);
   __sdag_pup(p);

    //p | u;
    //p | u2;
    //p | u3;
    p | imsg;
}
    
Advection::~Advection(){
    ckout << "In Destructor" << endl;
    delete [] u;
    delete [] u2;
    delete [] u3;

    delete [] x;
    delete [] y;
    delete [] left_edge;
    delete [] right_edge;
}

void Advection::begin_iteration(void) {
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
    
void Advection::process(int dir, int size, double gh[]){
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
    
void Advection::compute_and_iterate(){
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

void Advection::iterate() {
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
    
void Advection::requestNextFrame(liveVizRequestMsg *m){
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

#include "Advection.def.h"
