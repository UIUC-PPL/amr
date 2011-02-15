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
extern double xctr, yctr, radius;
extern double v;
extern double ap, an;
extern double tmax, t, dt, cfl;
extern int  ny;
#define index(i,j)  ((j)*(block_width+2) + i)

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

    mem_allocate(top_edge, block_width);
    mem_allocate(bottom_edge, block_width);
    
}

Advection::Advection(bool exists, bool isRefined, double xmin, double xmax, double ymin, double ymax){

    CBase_Advection();
    this->exists = true;
    this->isRefined = isRefined;
    
    for(int dir=UP; dir<=RIGHT; ++dir)
        nbr[dir] = thisIndex.getNeighbor(dir);
    if(thisIndex.nbits!=0)
        parent = thisIndex.getParent();
    
    /*Lets set the nieghbors status if this is the root node */
    if(thisIndex.nbits==0)
        for(int i=0; i<NUM_NEIGHBORS; i++){
            nbr_exists[i]=true;
            nbr_isRefined[i]=false;
        }

    hasReceived = *new set<int>();

    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;

    dx = (xmax - xmin)/double(nx);
    dy = (ymax - ymin)/double(ny);

    advection();
}
    
void Advection::printState(){
    QuadIndex qindex = thisIndex;
    ckout << "Printing Status of Node " << qindex.getIndexString() << endl;
    ckout << "exists: " << exists << endl;
    ckout << "isRefined: " << isRefined << endl;
    ckout << "nbr_exists: ";
    for(int j=0; j<NUM_NEIGHBORS; j++)
        ckout << nbr_exists[j] << ", ";
    ckout << endl << "nbr_isRefined: ";
    for(int j=0; j<NUM_NEIGHBORS; j++)
        ckout << nbr_isRefined[j] << ", ";
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
#if 0
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

    p|exists;
    p|isRefined;
    for(int i=0; i<NUM_NEIGHBORS; i++){
      p|nbr_exists[i];
      p|nbr_isRefined[i];
      p|nbr_dataSent[i];
    }
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
    ckout << "In Destructor" << endl;
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

void Advection::begin_iteration(void) {
    //ckout << "String: " << thisIndex.getIndexString() << endl;
    for(int i=0; i<NUM_NEIGHBORS; i++)
        nbr_dataSent[i]=false;
    
    hasReceived.clear();

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
    /* If neighbor exists and is not refined just send the data
       If neighbor does not exit, extrapolate the data and send to the correposnding neighbor of the parent
       If neighbor exists and is refined then do nothing, wait for the data from the neighbors children to arrive, 
            iterpolate the data and send the data */

    //ckout << "Left Neighbor of " << thisIndex.getIndexString() << " is " << nbr[LEFT].getIndexString() << endl;
    //send my left edge
    if(nbr_exists[LEFT] && !nbr_isRefined[LEFT]){ckout << "LEFT" << endl;
        thisProxy(nbr[LEFT]).receiveGhosts(iterations, RIGHT, block_height, left_edge);
    }
    else if(!nbr_exists[LEFT]){
        //extrapolate the data
        QuadIndex receiver = thisIndex.getParent().getNeighbor(LEFT);
        for(int j=1; j<=block_height; j+=2){
            left_edge[j/2] = (u[index(1,j)] + u[index(2,j)] + u[index(1,j+1)] +u[index(2,j+1)])/4;
        }
        thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), LEFT), block_height/2, left_edge);
    }
    
    //ckout << "Right Neighbor of " << thisIndex.getIndexString() << " is " << nbr[RIGHT].getIndexString() << endl;
    //send my right edge
    if(nbr_exists[RIGHT] && !nbr_isRefined[RIGHT]){ckout << "RIGHT" << endl;
        thisProxy(nbr[RIGHT]).receiveGhosts(iterations, LEFT, block_height, right_edge);
    }
    else if(!nbr_exists[RIGHT]){
        QuadIndex receiver = thisIndex.getParent().getNeighbor(RIGHT);
        for(int j=1; j<=block_height; j+=2){
            right_edge[j/2] = (u[index(block_width-1,j)] + u[index(block_width,j)] + u[index(block_width-1,j+1)] +u[index(block_width, j+1)])/4;
        }
        thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), RIGHT), block_height/2, left_edge);
    }

    //ckout << "Top Neighbor of " << thisIndex.getIndexString() << " is " << nbr[UP].getIndexString() << endl;
    //send my top edge
    if(nbr_exists[UP] && !nbr_isRefined[UP]){ckout << "UP: " << nbr[UP].getIndexString() <<endl;
        thisProxy(nbr[UP]).receiveGhosts(iterations, DOWN, block_width, &u[index(1,1)]);
    }
    else if(!nbr_exists[UP]){
        QuadIndex receiver = thisIndex.getParent().getNeighbor(UP);
        for(int i=1; i<=block_width; i+=2){
            top_edge[i/2] = (u[index(i,1)] + u[index(i,2)] + u[index(i+1,1)] + u[index(i+1,2)])/4;
        }
        thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), UP), block_width/2, top_edge);
    }
    
    //ckout << "Bottom Neighbor of " << thisIndex.getIndexString() << " is " << nbr[DOWN].getIndexString() << endl;
    //send my bottom edge
    if(nbr_exists[DOWN] && !nbr_isRefined[DOWN]){ckout << "DOWN" << endl;
        thisProxy(nbr[DOWN]).receiveGhosts(iterations, UP, block_width, &u[index(1, block_height)]);
    }
    else if(!nbr_exists[DOWN]){
        QuadIndex receiver = thisIndex.getParent().getNeighbor(DOWN);
        for(int i=1; i<=block_width; i+=2){
            bottom_edge[i/2] = (u[index(i,block_height-1)] + u[index(i,block_height)] + u[index(i+1,block_height-1)] + u[index(i+1,block_height)])/4;
        }
        thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), DOWN), block_width/2, bottom_edge);
    }

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
    ckout << "process called on " << dir << endl;
    switch(dir){
        case LEFT:
            imsg++;
            ckout << "Received From Left" << endl;
            hasReceived.insert(LEFT);
            for(int i=0; i<size; i++)
                u[index(0,i+1)] = gh[i];
                //u[i+1][0] = gh[i];
            break;

        case RIGHT:
            imsg++;
            ckout << "Received From RIGHT" << endl;
            hasReceived.insert(RIGHT);
            for(int i=0; i<size; i++)
                u[index(block_width+1,i+1)] = gh[i];
                //u[i+1][block_width+1]=gh[i];
            break;

        case UP:
            imsg++;
            ckout << "Received From UP" << endl;
            hasReceived.insert(UP);
            for(int i=0; i<size; i++)
                u[index(i+1,0)] = gh[i];
                //u[block_height+1][i+1]=gh[i];
            break;

        case DOWN:
            imsg++;
            ckout << "Received From Down" << endl;
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
    //check if data can be sent to any of the refined neighbors
    //If the neighbors are at the same level or do not exist at all 
    //data will be sent as a part of the current iteration and need 
    //not be sent here
    for(int i=0; i<NUM_NEIGHBORS; i++)
        if(nbr_isRefined[i] && !nbr_dataSent[i]){
            if(i==RIGHT){
                if(hasReceived.find(RIGHT_UP)!=hasReceived.end() &&
                    hasReceived.find(RIGHT_DOWN)!=hasReceived.end() && 
                    hasReceived.find(UP_RIGHT)!=hasReceived.end() &&
                    hasReceived.find(DOWN_RIGHT)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else if(i==LEFT){
                if(hasReceived.find(LEFT_UP)!=hasReceived.end() &&
                    hasReceived.find(LEFT_DOWN)!=hasReceived.end() &&
                    hasReceived.find(UP_LEFT)!=hasReceived.end() &&
                    hasReceived.find(DOWN_LEFT)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else if(i==UP){
                if(hasReceived.find(UP_LEFT)!=hasReceived.end() &&
                    hasReceived.find(UP_RIGHT)!=hasReceived.end() &&
                    hasReceived.find(LEFT_UP)!= hasReceived.end() &&
                    hasReceived.find(RIGHT_UP)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else{ // if i==DOWN
                if(hasReceived.find(DOWN_LEFT)!=hasReceived.end() &&
                    hasReceived.find(DOWN_RIGHT)!=hasReceived.end() &&
                    hasReceived.find(LEFT_DOWN) != hasReceived.end() &&
                    hasReceived.find(RIGHT_DOWN)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
        }
}
    
void Advection::interpolateAndSend(int NBR){
    double sx, sy;
    if(NBR==RIGHT){
        for(int i=0; i<block_height; i++){
            sx = (u[index(block_width-1, i+1)] - u[index(block_width+1, i+1)])/(2*dx);
            sy = (u[index(block_width, i)] - u[index(block_width, i+2)])/(2*dy);

            right_edge[wrap(2*i, block_height)] = u[index(block_width, i+1)] + sx*(dx/4) - sy*(dy/4);
            right_edge[wrap(2*i+1, block_height)] = u[index(block_width, i+1)] + sx*(dx/4) + sy*(dy/4);
            if(i==block_height/2 -1 ){//send the data to RIGHT_UP Neighbor
                QuadIndex receiver = nbr[RIGHT].getChild(map_child(LEFT_UP));//LEFT_UP child of the nieghbor
                thisProxy(receiver).receiveGhosts(iterations, LEFT, block_height, right_edge);
            }
            else if(i==block_height-1){// send the data to the RIGHT_DOWN Neighbor
                QuadIndex receiver = nbr[RIGHT].getChild(map_child(LEFT_DOWN));//LEFT_DOWN child of the neighbor
                thisProxy.receiveGhosts(iterations, LEFT, block_height, right_edge);
            }
        }
    }
    else if(NBR==LEFT){
        for(int i=0; i<block_height; i++){
            sx = (u[index(0, i+1)] - u[index(2,i+1)])/(2*dx);
            sy = (u[index(1, i)] - u[index(1, i+2)])/(2*dy);

            left_edge[wrap(2*i, block_height)] = u[index(1,i+1)] - sx*(dx/4) -sy*(dy/4);
            left_edge[wrap(2*i+1, block_height)] = u[index(1,i+1)] - sx*(dx/4) +sy*(dy/4);
        
            if(i==block_height/2 -1 ){//send the data to LEFT_UP Neighbor
                QuadIndex receiver = nbr[LEFT].getChild(map_child(RIGHT_UP));//LEFT_UP child of the nieghbor
                thisProxy(receiver).receiveGhosts(iterations, RIGHT, block_height, left_edge);
            }
            else if(i==block_height-1){// send the data to the LEFT_DOWN Neighbor
                QuadIndex receiver = nbr[LEFT].getChild(map_child(RIGHT_DOWN));//LEFT_DOWN child of the neighbor
                thisProxy.receiveGhosts(iterations, RIGHT, block_height, left_edge);
            }
        }
    }
    else if(NBR==UP){
        for(int i=0; i<block_width; i++){
            sx = (u[index(i,1)]-u[index(i+2,1)])/(2*dx);
            sy = (u[index(i+1,0)]-u[index(i+1,2)])/(2*dy);

            top_edge[wrap(2*i, block_width)] = u[index(i+1, 1)] -sx*(dx/4) - sy*(dy/4);
            top_edge[wrap(2*i+1, block_width)] = u[index(i+1, 1)] + sx*(dx/4) -sy*(dy/4);

            if(i==block_width/2-1){// send the data to UP_LEFT Neighbor
                QuadIndex receiver = nbr[UP].getChild(map_child(DOWN_LEFT));
                thisProxy(receiver).receiveGhosts(iterations, DOWN, block_width, top_edge);
            }
            else if(i==block_width-1){//send the data to the UP_RIGHT Neighbor
                QuadIndex receiver = nbr[UP].getChild(map_child(DOWN_RIGHT));
                thisProxy(receiver).receiveGhosts(iterations, DOWN, block_width, top_edge);
            }
        }
    }
    else{// if NBR==DOWN
        for(int i=0; i<block_width; i++){
            sx = (u[index(i,block_height)]-u[index(i+2, block_height)])/(2*dx);
            sy = (u[index(i+1, block_height-1)] - u[index(i+1, block_height+1)])/(2*dy);

            bottom_edge[wrap(2*i, block_width)] = u[index(i+1,block_height)] - sx*(dx/4) + sy*(dy/4);
            bottom_edge[wrap(2*i+1, block_width)] = u[index(i+1, block_height)] + sx*(dx/4) + sy*(dy/4);

            if(i==block_width/2-1){// send the data to UP_LEFT Neighbor
                QuadIndex receiver = nbr[DOWN].getChild(map_child(UP_LEFT));
                thisProxy(receiver).receiveGhosts(iterations, UP, block_width, bottom_edge);
            }
            else if(i==block_width-1){//send the data to the UP_RIGHT Neighbor
                QuadIndex receiver = nbr[DOWN].getChild(map_child(UP_RIGHT));
                thisProxy(receiver).receiveGhosts(iterations, UP, block_width, bottom_edge);
            }
        }
    }
};

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
    
#if 1
    for(int i=1; i<=block_height; i++){
        for(int j=1; j<=block_width; j++)
            ckout << u[index(j,i)] << "\t";
            ckout << endl;
    }
    ckout << endl;
    ckout << "After First Iteration" << endl;
#endif
    iterate();
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
    if (isRefined){
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

void Advection::interpolate(double *u, double *refined_u, int xstart, int xend, int ystart, int yend){

    for(int i=xstart; i<=xend; i++)
        for(int j=ystart; j<=yend; j++){
            sx = (u[index(i-1,j)]-u[index(i+1,j)])/(2*dx);
            sy = (u[index(i,j-1)]-u[index(i,j+1)])/(2*dy);

            refined_u[index(2*(i-1), 2*(j-1))] = u[index(i,j)] - sx*(dx/4) - sy*(dy/4);
            refined_u[index(2*(i-1)+1, 2*(j-1))] = u[index(i,j)] + sx*(dx/4) - sy*(dy/4);
            refined_u[index(2*(i-1), 2*(j-1)+1)] = u[index(i,j)] - sx*(dx/4) + sy*(dy/4);
            refined_u[index(2*(i-1)+1, 2*(j-1)+1)] = u[index(i,j)] + sx*(dx/4) + sy*(dy/4);
        }
}

void Advection::refine(){
    //Spawn The four children and give them the data
    //Assuming we already have the new boundary data
    
    //Interpolate the data and give it to the children when they are initialized
    // boundaries of the children will have to be sent by the neighbor
    double sx, sy, s;
    double *refined_u = new double[(block_width+2)*(block_height+2)];

    interpolate(u, refined_u, 1, block_width/2, 1, block_height/2);
    //initialize the child
    thisProxy(thisIndex.getChild("01")).insert(xmin, xmax, ymin, ymax, refined_u);

    interpolate(u, refined_u, block_width/2+1, block_width, 1, block_height/2);
    //initialize the child

    interpolate(u, refined_u, 1, block_width/2, block_height/2+1, block_height);
    //init the child

    interpolate(u, refined_u, block_width/2+1, block_width, block_height/2+1, block_height);
    //init the child
}

void Advection::Advection(double, double, double, double, double*){

}

#include "Advection.def.h"
