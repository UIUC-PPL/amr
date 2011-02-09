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
extern CProxy_Main mainProxy;

extern int array_height;
extern int array_width;

extern int num_chare_rows;
extern int num_chare_cols;

extern int block_height;
extern int block_width;

#define wrap_x(a) (((a)+num_chare_rows)%num_chare_rows)
#define wrap_y(a) (((a)+num_chare_cols)%num_chare_cols)


extern CkArrayID a;

extern int nframe;
extern double xmin, xmax, ymin, ymax;
extern double xctr, yctr, radius;
extern int nx, ny;
extern double dx, dy, v;
extern double ap, an;
extern double tmax, t, dt, cfl;


#define index(i,j)  ((j)*(block_width+2) + i)


void Advection::mem_allocate(double* &p, int size){
    p = new double[size];
}


void Advection::mem_allocate_all(){

    mem_allocate(left_edge, block_height);
    mem_allocate(right_edge, block_height);

    //if(!isdummy){
        mem_allocate(u, (block_width+2)*(block_height+2));
        mem_allocate(u2, (block_width+2)*(block_height+2));
        mem_allocate(u3, (block_width+2)*(block_height+2));

        mem_allocate(x, block_width+2);
        mem_allocate(y, block_height+2);
    //}
   // else{
        mem_allocate(top_edge, block_width);
        mem_allocate(bottom_edge, block_width);
    //}
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

    advection();
}
    
    
void Advection::refine(){	
    // should always be called on a node not having Real Chilren
    //ckout << "Refine called on " << thisIndex.getIndexString() << endl;

    if(hasRealChildren)//somehow if it has already (refined due to some other call) then do nothing
        return;

    if(hasDummyChildren){               
        hasDummyChildren = false;
        hasRealChildren = true;

        for(int i=0; i<NUM_CHILDREN; i++){
            thisProxy(thisIndex.getChild(i)).setReal();
        }
    }
    else if(!hasRealChildren){
        // Create new real chilren
        hasDummyChildren = false;
        hasRealChildren = true;

        for(int i=0; i<4; i++){
            thisProxy(thisIndex.getChild(i)).insert(false, false, false); 
        }
        thisProxy.doneInserting();
    }
    
    for(int i=0; i<NUM_NEIGHBORS; i++){
        //if(thisIndex!=nbr[i])
        thisProxy(nbr[i]).inform_nbr_of_refinement(SENDER_DIR[i]);
    }
}

void Advection::inform_nbr_of_refinement(int inbr){
    /*Update Neighbor Infomration*/
    //ckout << inbr << " has refined" << endl;
    nbr_isdummy[inbr]=false;
    nbr_hasRealChildren[inbr]=true;
    nbr_hasDummyChildren[inbr]=false;
    
    /* Take Action To Maintain Mesh Structure*/
    if(isdummy){
        /*Two things need to be done here:
         * 1. Ask the parent to refine itself so that myself and my siblings become real.
         * 2. Create my dummy children to be available for communication with the neighbor who has decided to refine itself
         */
        thisProxy(parent).refine();
    }    
    else if(!hasDummyChildren && !hasRealChildren){
        hasDummyChildren = true;/* Will Now create dummy children*/
        hasRealChildren = false;
    }
    else if(hasRealChildren || hasDummyChildren)
        return;

    //create dummy children  
    for(int i=0; i<4; i++)
        thisProxy(thisIndex.getChild(i)).insert(true, false, false);      

    thisProxy.doneInserting();
    /* Inform Neighbors that I now have Dummy Children*/
    for(int i=0; i<NUM_NEIGHBORS; i++){
        thisProxy(nbr[i]).inform_nbr_hasDummyChildren(SENDER_DIR[i]);
    }
    /* Inform Parent Also that I now have Dummy Children */
    thisProxy(parent).inform_child_hasDummyChildren(thisIndex.getChildNum());
}
 
void Advection::derefine(){
    /* If Any of My Children Have Dummy Children then I cannot derefine
     * because that means the dummy grandchildren are serving some of their
     * real neighbors and hence cannot be removed*/
    for(int i=0; i<NUM_CHILDREN; i++)
        if(child_hasDummyChildren[i]==true){
            ckout << "Cannot derefine MySelf" << endl;
            return;
        }
    /* Now I can Derefine MySelf */

    /* Check If Any of My Neighbors have Real Children Then I have to make 
     * my children Dummy */
    int i=0;
    for(i=0; i<NUM_NEIGHBORS; i++){
        if(nbr_hasRealChildren[i]){
            /* Make my Children Dummy*/
            for(int j=0; j<NUM_CHILDREN;j++)
                thisProxy(thisIndex.getChild(j)).setDummy();
            /* Inform Neighbors that I now have Dummy  Children*/
            for(int j=0; j<NUM_NEIGHBORS; j++)
                thisProxy(nbr[j]).inform_nbr_hasDummyChildren(SENDER_DIR[j]);
            return;
        }
    }
    if (i==NUM_NEIGHBORS){
        /* i.e. when none of the neighbors have real children
         * In that case, delete the children
         */
        destroyChildren();
        /* Inform Neighbors of What You Did to Your Children*/
        for(int j=0; j<NUM_NEIGHBORS; j++)
            thisProxy(nbr[j]).inform_nbr_hasNoChildren(SENDER_DIR[j]);
    }
}

void Advection::destroyChildren(){
    for(int j=0; j<NUM_CHILDREN; j++){
        // free the Memory Held by the Children
        thisProxy(thisIndex.getChild(j)).free_memory();
        //destoy the array element
        thisProxy(thisIndex.getChild(j)).ckDestroy();
    }
}

void Advection::inform_nbr_hasNoChildren(int inbr){
    /* This entry method will be called when the neighbor has derefined 
     * and has destroyed his children, the neighbor can derefine in such a 
     * manneronly if I had no Real Children. That means I can have Dummy Children.
     * Now I should check that if I have Dummy Children and there are no neighbors 
     * with Real Children, then I should also destry my dummy children*/
     
    nbr_hasDummyChildren[inbr]=false;
    nbr_hasRealChildren[inbr]=false; 
    int i;
    if(hasDummyChildren){/* do not destroy if it has real children
                          * because may be those children have been created
                          * as a result of refinement that occured simultaneously with
                          * neighbors derefinement */
        for(i=0; i<NUM_NEIGHBORS; i++){
            if(nbr_hasRealChildren[i]==true)
                /* If any of the other neighbors has real children
                 * then donot destroy your dummy children*/                                            
                break;
        }
        if(i==NUM_NEIGHBORS)
            destroyChildren();
    }
}

void Advection::inform_nbr_hasDummyChildren(int inbr){
    nbr_hasDummyChildren[inbr]=true;
    /* Check if all the neighbors have dummy children
     * And if That is the case destroy all Children */
    int i;
    if(hasDummyChildren){/* If I have Dummy Children then I can destroy them 
                          * if none of my neighbors has real children */  
        for(i=0; i<NUM_NEIGHBORS; i++){
            if(nbr_hasRealChildren[i]==true)
                break;
        }
        if(i==NUM_NEIGHBORS)
            destroyChildren();
    }
}

void Advection::printState(){
    QuadIndex qindex = thisIndex;
    ckout << "Printing Status of Node " << qindex.getIndexString() << endl;
    ckout << "isdummy: " << isdummy << endl;
    ckout << "hasRealChildren: " << hasRealChildren << endl;
    ckout << "hasDummyChldren: " << hasDummyChildren << endl;
    ckout << "nbr_isdummy: ";
    for(int j=0; j<NUM_NEIGHBORS; j++)
        ckout << nbr_isdummy[j] << ", ";
    ckout << endl << "nbr_hasRealChildren: ";
    for(int j=0; j<NUM_NEIGHBORS; j++)
        ckout << nbr_hasRealChildren[j] << ", ";
    ckout << endl << "nbr_hasDummyChildren: ";
    for(int j=0; j<NUM_NEIGHBORS; j++)
        ckout << nbr_hasDummyChildren[j] << ", ";
    ckout << endl << "child_hasDummyChildren: ";
    for(int j=0; j<NUM_CHILDREN; j++)
        ckout << child_hasDummyChildren[j] << ", ";
    ckout << endl;
}

void Advection::advection(){
    __sdag_init();
    usesAtSync = CmiTrue;
    //ckout << "In Advection: " << endl;
    myt = t;
    mem_allocate_all();
    iterations=0;

    thisIndex.getCoordinates(xc, yc);

    for(int i=0; i<block_width+2; i++){
        x[i] = xmin + double(xc*block_width+i)*dx - 0.5*dx;
        //ckout << x[i] << endl;
    }
    //ckout << endl;

    for(int i=0; i<block_height+2; i++){
        y[i] = ymin + double(yc*block_height+i)*dy - 0.5*dy;
       // ckout << y[i] << endl;
    }
    //ckout << endl;

    double rsq;
    
    //ckout << "In Adfvection2" << endl;
    //ckout << xctr << ", " << yctr << endl;
    for(int i=0; i<block_height+2; i++){
        for(int j=0; j<block_width+2; j++){
            rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[j]-yctr);
            if(rsq <= radius*radius)
                u[index(i,j)] = 2;
            else u[index(i,j)] = 1;
        }
    }
    //ckout << "In Advection 3: " <<u[1][1]<< endl;
#if 1
    for(int i=0; i<block_height; i++){
        for(int j=0; j<block_width; j++)
            ckout << u[index(j+1,i+1)] << "\t";
        ckout << endl;
    }
    ckout << endl;
#endif
}

    //added for array migration - see how 2D arrays can be packed
void Advection::pup(PUP::er &p){
    ckout << "In PUP" << endl;
    CBase_Advection::pup(p);
   __sdag_pup(p);

    p|isdummy;
    p|hasRealChildren;
    p|hasDummyChildren;
    for(int i=0; i<NUM_NEIGHBORS; i++)
      p|nbr[i];
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
    }

    for (int i=0; i<block_height+2; i++){
        p|y[i];
    }
    for (int i=0; i<block_height+2; i++){
        p|left_edge[i];
        p|right_edge[i];
    }

    p|iterations;
    p|up;
    p|un;
    p|myt;
    p|mydt;
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

    delete [] top_edge;
    delete [] right_edge;
}

void Advection::setReal(){
    isdummy = false;
    manage_memory_DummyToReal();
}

void Advection::setDummy(){
    isdummy = true;
    manage_memory_RealToDummy();
}

void Advection::manage_memory_RealToDummy(){
    delete [] u;
    delete [] u2;
    delete [] u3;
    
    delete [] x;
    delete [] y;

    mem_allocate(top_edge, block_width+2);
    mem_allocate(bottom_edge, block_width+2);
}


void Advection::manage_memory_DummyToReal(){
    delete [] top_edge;
    delete [] bottom_edge;
    
    mem_allocate(u, (block_height+2)*(block_width+2));
    mem_allocate(u2, (block_height+2)*(block_width+2));
    mem_allocate(u3, (block_height+2)*(block_width+2));

    mem_allocate(x, block_width+2);
    mem_allocate(y, block_height+2);
}

void Advection::begin_iteration(void) {
    //ckout << "String: " << thisIndex.getIndexString() << endl;

    for(int i=1; i<=block_width; i++)
        for(int j=1; j<=block_height; j++)
            u2[index(i,j)] = u[index(i,j)];

    for(int i=1; i<=block_width; i++)
        for(int j=1; j<=block_height; j++)
            u3[index(i,j)] = u[index(i,j)];

    iterations++;

    for(int j=1; j<=block_height; j++){
        left_edge[j-1] = u[index(1,j)];
        right_edge[j-1] = u[index(block_width,j)];
    }
        ///*
    //ckout << "Left Neighbor of " << thisIndex.getIndexString() << " is " << nbr[LEFT].getIndexString() << endl;
    //send my left edge
    thisProxy(nbr[LEFT]).receiveGhosts(iterations, RIGHT, block_height, left_edge);

    //ckout << "Right Neighbor of " << thisIndex.getIndexString() << " is " << nbr[RIGHT].getIndexString() << endl;
    //send my right edge
    thisProxy(nbr[RIGHT]).receiveGhosts(iterations, LEFT, block_height, right_edge);


    //ckout << "Top Neighbor of " << thisIndex.getIndexString() << " is " << nbr[UP].getIndexString() << endl;
    //send my top edge
    thisProxy(nbr[UP]).receiveGhosts(iterations, DOWN, block_width, &u[index(1,1)]);

    //ckout << "Bottom Neighbor of " << thisIndex.getIndexString() << " is " << nbr[DOWN].getIndexString() << endl;
    //send my bottom edge
    thisProxy(nbr[DOWN]).receiveGhosts(iterations, UP, block_width, &u[index(1, block_height)]);

}
template<class T>
void print_Array(T* array, int size, int row_size){
    for(int i=0; i<size; i++){
        if(i%row_size==0)
            ckout << endl;
        ckout << array[i] << '\t';
    }
}

void Advection::process(int iter, int dir, int size, double gh[]){
//printf("[%d] process %d %d\n", thisIndex, iter, dir);
    switch(dir){
        case LEFT:
            for(int i=0; i<size; i++)
                u[index(0,i+1)] = gh[i];
                //u[i+1][0] = gh[i];
        break;

        case RIGHT:
            for(int i=0; i<size; i++)
                u[index(block_width+1,i+1)] = gh[i];
                //u[i+1][block_width+1]=gh[i];
        break;

        case UP:
            for(int i=0; i<size; i++)
                u[index(i+1,0)] = gh[i];
                //u[block_height+1][i+1]=gh[i];
            break;

         case DOWN:
            for(int i=0; i<size; i++)
                u[index(i+1,block_height+1)] = gh[i];
                //u[0][i+1]=gh[i];
            break;
        default:
            CkAbort("ERROR\n");
    }
}
    
void Advection::compute_and_iterate(){
    //ckout << "dt: " << dt << " ap:" << ap << " an:" << an << endl;
    //ckout << iterations << endl;
    for(int i=1; i<=block_width; i++){
        for(int j=1; j<=block_height; j++){
            up = (u[index(i+1,j)] - u[index(i,j)])/dx;
            un = (u[index(i,j)]-u[index(i-1,j)])/dx;

            u2[index(i,j)] = u[index(i,j)] - dt* (ap*un + an*up);
            
            //ckout << u2[index(i,j)] << endl;
        }
    }

    for(int i=1; i<=block_width; i++)
        for(int j=1; j<=block_height; j++){
            up = (u[index(i,j+1)] - u[index(i,j)])/dy;
            un = (u[index(i,j)] - u[index(i,j-1)])/dy;

            u3[index(i,j)] = u[index(i,j)] - dt*(ap*un + an*up);

            //ckout << "u3:" << u3[index(i,j)] << endl;
        }

    for(int j=1; j<=block_height; j++)
        for(int i=1; i<=block_width; i++)
           u[index(i,j)] = 0.5*(u2[index(i,j)] + u3[index(i,j)]);
    
    //Check Whether it is the time for Refinement
    /*if(checkForRefinement()){
        refine();
        return;
    }
    else if(checkForDerefinement()){
        derefine();
        return;
    }*/
    

#if 1
    for(int i=1; i<=block_height; i++){
        for(int j=1; j<=block_width; j++)
            ckout << u[index(j,i)] << "\t";
            ckout << endl;
    }
    ckout << endl;
    ckout << "After Iteration " << iterations <<  " on " << thisIndex.getIndexString() << endl;
#endif
    iterate();
}

bool Advection::checkForRefinement(){
    if(rand()%10==0)
        return true;
    return false;
}

bool Advection::checkForDerefinement(){
    if(rand()%10==0)
        return true;
    return false;
}

void Advection::iterate() {
         myt = myt+mydt;
         if(myt < tmax){
             mydt = min(dx,dy)/v * cfl;
             if ((myt + mydt) >= tmax )
                 mydt = tmax - myt;
             doStep();    
         }
         else {
             CkPrintf("Contribute\n");
             contribute();
         }
}

void Advection::requestNextFrame(liveVizRequestMsg *m){
    if (isdummy || hasRealChildren){
        liveVizDeposit(m,0,0,0,0,NULL, this);
    }

    int sy = yc;//thisIndex%num_chare_cols;
    int sx = xc;//thisIndex/num_chare_rows;
    int w= block_width, h = block_height;

    unsigned char *intensity = new unsigned char[3*w*h];
    for(int j=0; j<h; j++){
       for(int i=0; i<w; i++){
           // ckout << u[i+1][j+1] << endl;
            intensity[3*(i+j*w)+ 0] = 255;// red component
            intensity[3*(i+j*w)+ 1] = 255 - (u[index(i+1,j+1)]-1)*255;// BLUE component
            intensity[3*(i+j*w)+ 2] = 255 - (u[index(i+1,j+1)]-1)*255;// green component
        }
    }
    liveVizDeposit(m, sy*w,sx*h, w,h, intensity, this);
    delete[] intensity;
}
#include "Advection.def.h"
