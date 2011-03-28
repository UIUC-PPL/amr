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
extern int nx, ny;
extern map<string, DIR> nbrDirectionMap;
extern map<DIR, DIR> reverse_dir_map;

#define index(i,j)  (int)((j)*(block_width+2) + i)


InitRefineMsg::InitRefineMsg(double dx, double dy, double myt, double mydt, int iterations, double *refined_u, bool nbr_exists[NUM_NEIGHBORS], bool nbr_isRefined[NUM_NEIGHBORS], DECISION nbr_decision[3*NUM_NEIGHBORS]){
    this->dx = dx;
    this->dy = dy;
    this->myt = myt;
    this->mydt = mydt;
    this->iterations = iterations;
    memcpy(this->refined_u, refined_u, sizeof(double)*block_height*block_width);
    memcpy(this->parent_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
    memcpy(this->parent_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
    memcpy(this->parent_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);
}

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

Advection::Advection(double xmin, double xmax, double ymin, double ymax){
//Constructor for the Initial Grid Zones
    CBase_Advection();
    this->exists = true;
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
    ckout << "xctr: " << xctr << ", yctr: " << yctr << ", radius: " << radius << endl;
    for(int i=0; i<block_height+2; i++){
        for(int j=0; j<block_width+2; j++){
            rsq = (x[i] - xctr)*(x[i]-xctr) + (y[j] - yctr)*(y[j]-yctr);
            if(rsq <= radius*radius)
                u[index(i,j)] = 2;
            else u[index(i,j)] = 1;
        }
    }
#if 1
    for(int i=0; i<block_height; i++){
        for(int j=0; j<block_width; j++)
            ckout << u[index(j+1,i+1)] << "\t";
        ckout << endl;
    }
    ckout << endl;
    //CkExit();
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

void Advection::sendGhost(int dir){
    if(nbr_exists[dir] && !nbr_isRefined[dir]){
        ckout << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Neighbor " << nbr[dir].getIndexString() << ", iteration " << iterations << endl;
        if(dir==LEFT)
            thisProxy(nbr[dir]).receiveGhosts(iterations, SENDER_DIR[dir], block_height, left_edge);
        else if(dir==RIGHT)
            thisProxy(nbr[dir]).receiveGhosts(iterations, SENDER_DIR[dir], block_height, right_edge);
        else if (dir==UP)
            thisProxy(nbr[dir]).receiveGhosts(iterations, SENDER_DIR[dir],  block_width, &u[index(1,1)]);
        else if(dir==DOWN)
              thisProxy(nbr[dir]).receiveGhosts(iterations, SENDER_DIR[dir],  block_width, &u[index(1, block_height)]);
    }
    else if(!nbr_exists[dir]){
        QuadIndex receiver = thisIndex.getParent().getNeighbor(dir);
        ckout << thisIndex.getIndexString() << " sending Ghost in Dir " << dir << " to Uncle: " << receiver.getIndexString() << ", iteration " << iterations << endl;
        if(dir==LEFT){
            for(int j=1; j<=block_height; j+=2){
                left_edge[j/2] = (u[index(1,j)] + u[index(2,j)] + u[index(1,j+1)] +u[index(2,j+1)])/4;
            }
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), LEFT), block_height/2, left_edge);
        }
        else if(dir==RIGHT){
            for(int j=1; j<=block_height; j+=2){
                right_edge[j/2] = (u[index(block_width-1,j)] + u[index(block_width,j)] + u[index(block_width-1,j+1)] +u[index(block_width, j+1)])/4;
            }
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), RIGHT), block_height/2, left_edge);
        }
        else if(dir==UP){
            for(int i=1; i<=block_width; i+=2){
                top_edge[i/2] = (u[index(i,1)] + u[index(i,2)] + u[index(i+1,1)] + u[index(i+1,2)])/4;
            }
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), UP), block_width/2, top_edge);
        }
        else if(dir==DOWN){
            for(int i=1; i<=block_width; i+=2){
                bottom_edge[i/2] = (u[index(i,block_height-1)] + u[index(i,block_height)] + u[index(i+1,block_height-1)] + u[index(i+1,block_height)])/4;
            }
            thisProxy(receiver).receiveGhosts(iterations, map_nbr(thisIndex.getChildNum(), DOWN), block_width/2, bottom_edge);
        }
    }
}

void Advection::begin_iteration(void) {
    //ckout << "String: " << thisIndex.getIndexString() << endl;
    ckout << "Begin Iteration " << iterations << " on " << thisIndex.getIndexString() << endl;
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
    //ckout << "Seding Ghosts " << thisIndex.getIndexString() << endl;
    for(int i=0; i<NUM_NEIGHBORS; i++){
        sendGhost(i);
    }
    ckout << "Done Sending Ghosts " << thisIndex.getIndexString() << endl;
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
    ckout << thisIndex.getIndexString() << " received data for direction " << dir << ", iteration " << iter << endl;
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
}

void Advection::sendReadyData(){
    //check if data can be sent to any of the refined neighbors
    //If the neighbors are at the same level or do not exist at all 
    //data will be sent in begin_iteration function and need 
    //not be sent here
    for(int i=0; i<NUM_NEIGHBORS; i++)
        if(nbr_isRefined[i] && !nbr_dataSent[i]){
            if(i==RIGHT){
                if(hasReceived.find(RIGHT_UP)!=hasReceived.end() &&
                    hasReceived.find(RIGHT_DOWN)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else if(i==LEFT){
                if(hasReceived.find(LEFT_UP)!=hasReceived.end() &&
                    hasReceived.find(LEFT_DOWN)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else if(i==UP){
                if(hasReceived.find(UP_LEFT)!=hasReceived.end() &&
                    hasReceived.find(UP_RIGHT)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
            else{ // if i==DOWN
                if(hasReceived.find(DOWN_LEFT)!=hasReceived.end() &&
                    hasReceived.find(DOWN_RIGHT)!=hasReceived.end()){
                    interpolateAndSend(i);
                    nbr_dataSent[i]=true;
                }
            }
        }
        /*if(imsg==4)
            compute_and_iterate();*/
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
	     if(iterations%10==0){//time to check need for refinement/coarsening
	     	/*This computation phase can be tested for correctness by running 
	     	extreme cases - like everyone wants to refine, 
	     	everyone wants to derefine, nobody wants to do anything */
		hasInitiatedPhase1=false;
                ckout << "Entering Mesh Restructure Phase on " << thisIndex.getIndexString() << endl;
		doMeshRestructure();
	     }
	     else{
             	doStep();    
	     }
         }
         else {
             CkPrintf("Contribute\n");
             contribute();
         }
}

DECISION Advection::getGranularityDecision(){
    if(strcmp(thisIndex.getIndexString(),"00")==0)// && iterations <15)
        return REFINE;
    else return STAY;
    /*if(thisIndex.getDepth()==min_depth)
        return STAY;*/
    

    /*if(rand()%3==0){
        return REFINE;
    }
    else{
        if(rand()%3==0 && thisIndex.getDepth()!=min_depth))
	    return DEREFINE;	
	else 
	    return STAY;
    }*/
}

void Advection::doMeshRestructure(){
    if(!hasInitiatedPhase1){
        hasInitiatedPhase1=true;

        if(!isRefined){//run this on leaf nodes
            /*Reset Old Data*/
            /*Phase1 Resetting*/
            for(int i=0; i<3*NUM_NEIGHBORS; i++)
                nbr_decision[i]=INV;

            decision = INV;

            for(int i=0; i<NUM_CHILDREN; i++)
                child_decision[i]=INV;

            for(int i=0; i<NUM_NEIGHBORS; i++)
                hasReceivedStatusFromParent[i]=false;

            parentHasAlreadyMadeDecision=false;
            hasReceivedParentDecision=false;

            /*Phase2 resetting*/
            hasInitiatedPhase2 = false;
            for(int i=0; i<NUM_NEIGHBORS; i++)
                nbr_dataSent[i]=false;

            hasReceived.empty();
            hasAllocatedMemory=false;

            /*Done Resetting Old Data*/
            if(!isRefined)
                decision = getGranularityDecision();

            //initiate Phase1 of the computation
            if(decision==REFINE||decision==STAY){
                communicatePhase1Msgs();
            }
            //Inform Parent About Start of the Restructure Phase 
            // Think of a better way to do this because right now as many as 
            // number of unrefined children are informing the parent about the
            // start of the restructure phase while only one of them need to do so
            if(parent!=thisIndex)
                thisProxy(parent).doMeshRestructure();
        }
	else if(isGrandParent && !parentHasAlreadyMadeDecision){
	    informParent(-1, INV);
	}

	//Wait for Quiescence
        
	CkStartQD(CkIndex_Advection::doPhase2(), &thishandle);
    }
}

/***** PHASE1 FUNCTIONS****/
void Advection::communicatePhase1Msgs(){

    if(decision==REFINE){
	//tell the neighbor if it exists, also tell to children of the neighbor 
	// if that neighbor is refined
	//In case the neighbor does not exist, just tell to the parent of the 
	//non-existing neighbor

	for(int i=0; i<NUM_NEIGHBORS; i++){
	    if(nbr_exists[i] && !nbr_isRefined[i])
	        thisProxy(nbr[i]).exchangePhase1Msg(SENDER_DIR[i]);//Since Phase1Msgs are only refinement messages
						      //just send your direction w.r.t. to the receiving neighbor
	    else if(nbr_exists[i] && nbr_isRefined[i]){
		    //Get Corresponding Children of the neighbor
		    QuadIndex q1, q2;
                    getChildren(nbr[i], SENDER_DIR[i], q1, q2);
                    thisProxy(q1).exchangePhase1Msg(SENDER_DIR[i]);
                    thisProxy(q2).exchangePhase1Msg(SENDER_DIR[i]);
	    }
	    else{//send to the parent of the non-existing neighbor
                thisProxy(nbr[i].getParent()).exchangePhase1Msg(map_nbr(nbr[i].getQuadI(), i));
	    }
	}
    }
    else{//No Need to do anything, just wait for 
    	 //neighbors to tell if they wish to derefine
    }

    //If my DECISION is to stay or to REFINE, tell the parent
    if((decision == REFINE || decision == STAY) && (parent!=thisIndex)){
    	thisProxy(parent).informParent(thisIndex.getChildNum(), decision);
    }
}

void Advection::informParent(int childNum, DECISION dec){//Run on Parent All of Whose Children are Leaf Nodes
    if(dec==REFINE){
        child_isRefined[childNum]=true;
        isGrandParent = true;
    }
    if(parentHasAlreadyMadeDecision==false){
        parentHasAlreadyMadeDecision=true;
	//tell rest of the children which are not refined
	for(int i=0 ;i<NUM_CHILDREN; i++){
	    if(i!=childNum && !child_isRefined[i])
	        thisProxy(thisIndex.getChild(i)).recvParentDecision();
	}
	//now tell the neighbors that I am Not Going to Derefine
	for(int i=0; i<NUM_NEIGHBORS; i++){//all the neighbors must be existing as I am a Parent
	    CkAssert(nbr_exists[i]);
	    thisProxy(nbr[i]).recvNeighborDecision(SENDER_DIR[i]);
	}
    }
}

void Advection::recvParentDecision(){
    hasReceivedParentDecision=true;
    if(decision==DEREFINE){
        decision=STAY;
    }
}

void Advection::recvNeighborDecision(DIR dir){
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
    else if(isRefined){/*then the parent already knows that it cannot derefine
                        however, tell the children that the neighbor is not going to derefine
                        */
        int c1, c2;
        thisIndex.getSiblingInDirection(dir, c1, c2);
        thisProxy(thisIndex.getChild(c1)).recvStatusUpdateFromParent(dir);
        thisProxy(thisIndex.getChild(c2)).recvStatusUpdateFromParent(dir);
    }
}

void Advection::recvStatusUpdateFromParent(int dir){
    hasReceivedStatusFromParent[dir]=true;
    if(nbr_decision[dir]!=REFINE){
        nbr_decision[dir]=STAY;
    }
}

void Advection::exchangePhase1Msg(int dir){//Phase1 Msgs are all Refine Messages
    //save the nbr's decision
    nbr_decision[dir]=REFINE;
    
    if(decision!=REFINE){
	//Now check if my Decision Changes Because of this Message
	if(dir==LEFT||dir==RIGHT||dir==UP|dir==DOWN);//don't do anything
	else if(dir==LEFT_UP||dir==LEFT_DOWN||dir==RIGHT_UP||dir==RIGHT_DOWN||
	    dir==UP_LEFT||dir==UP_RIGHT||dir==DOWN_LEFT||dir==DOWN_RIGHT){
	    decision = REFINE;
	    //tell all neighbors
	    communicatePhase1Msgs();
	}
    }
}


/**** PHASE2 FUNCTIONS ****/
void Advection::doPhase2(){
    hasInitiatedPhase1 = false; //reset this for the next MeshRestructure Phase for a Parent
    ckout << "Entering Phase2 on " << thisIndex.getIndexString() << endl;
    //Do Some Sanity Checks
    if(decision == DEREFINE){
        for(int i=0; i<NUM_NEIGHBORS; i++){
	    if(nbr_exists[i] && !nbr_isRefined[i] && nbr_decision[i]==REFINE){
	        ckout << "ERROR: doPhase2()::Sanity Check: My decision is Derefine But one of my\
		          Neighbor wants to refine" << endl;
		CkExit();
	    }
	    if(nbr_exists[i] && nbr_isRefined[i]){
	        if(nbr[i].getParent() == parent){//if same parent
		    ckout << "ERROR: doPhase()::Sanity Check: My decision is to derefine but one\
		    	      of my siblings has children" << endl;
                    CkExit();
		}
		else{//different parent
		   if(nbr_decision[i]!=STAY){
		       ckout << "ERROR: doPhase2::Sanity Check: MyDecision is to DEREFINE but one of\
		                 my refined neighbors hasn't yet said that it will derefine" << endl;
			CkExit();
		   }
		}
	    }
	}
    }

    //Update the decision of your neighbors
    if(decision == STAY){
        for(int i=0; i<NUM_NEIGHBORS; i++){
            if(!nbr_exists[i]){
                if(nbr[i].getParent().getDepth()==min_depth){
                    nbr_decision[i]=STAY;
                }
            }
        }
    }
    else if(decision==DEREFINE){// I am going to destroy myself
        for(int i=0; i<NUM_NEIGHBORS; i++){
            if(nbr_exists[i] && nbr_isRefined[i] && nbr_decision[i]!=STAY){
                //nbr_isRefined[i]=false;
                nbr_decision[i]=DEREFINE;
            }
            else if(nbr_exists[i] && !nbr_isRefined[i] && !hasReceivedStatusFromParent[i]){
                //nbr_exists[i]=false;
                nbr_decision[i]=DEREFINE;
            }
            else if(!nbr_exists[i] && nbr_decision[i]==REFINE){
                //nbr_exists[i]=true;
            }
        }
    }

    //Send Appropriate Data to the Neighbors based on their New Status
    if(isRefined==false){// can send data only if I am a leaf node
        for(int i=0; i<NUM_NEIGHBORS; i++){
            if(nbr_decision[i]==REFINE){
                if(i==LEFT){
                    for(int j=1; j<=block_height; j++)
                        left_edge[j-1] = u[index(1,j)];
                }
                else if(i==RIGHT){
                    for(int j=1; j<=block_height; j++)
                        right_edge[j-1] = u[index(block_width,j)];
                }
                sendGhost(i);
	    }
        }
    }

    if(isRefined && !isGrandParent && !parentHasAlreadyMadeDecision){//I am a parent(whose None of the Children Are Refined) and has to derefine
        //Get Data From the Children and extrapolate it
    }
    else if(decision==DEREFINE && !isRefined && !hasReceivedParentDecision){//send data to the parent
        size_t sz = ((block_height)*(block_width))/4;
        ChildDataMsg *msg = new (sz) ChildDataMsg();
    	//extrapolate the data
	for(int i=1; i<= block_width; i+=2){
	    for(int j=1; j<=block_height; j+=2){
	        int idx = index(i/2, j/2);
	        msg->child_u[idx] = u[index(i,j)];
		msg->child_u[idx] += u[index(i+1,j)];
		msg->child_u[idx] += u[index(i, j+1)];
		msg->child_u[idx] += u[index(i+1,j+1)];
		msg->child_u[idx] /= 4;
	    }
	}
	msg->childNum = thisIndex.getChildNum();
        memcpy(msg->child_nbr_exists, nbr_exists, sizeof(bool)*NUM_NEIGHBORS);
        memcpy(msg->child_nbr_isRefined, nbr_isRefined, sizeof(bool)*NUM_NEIGHBORS);
        memcpy(msg->child_nbr_decision, nbr_decision, sizeof(DECISION)*3*NUM_NEIGHBORS);

        thisProxy(parent).recvChildData(msg);
	//deallocate all your memory and destroy yourself
	thisProxy(thisIndex).ckDestroy();
        this->~Advection();
    }
    else if(decision==REFINE){
        getGhostsAndRefine();
    }
    
    //Update the Status of Your Neighbors
    if(decision == STAY){
        ckout << "Phase2: " << thisIndex.getIndexString() << " updating the Status of Neighbors" << endl;
       for(int i=0; i<NUM_NEIGHBORS; i++){
           if(!nbr_exists[i]){
	       if(nbr_decision[i]==DEREFINE){
	           ckout << "ERROR: doPhase(): Uncle Cannot DEREFINE while I want to STAY" << endl;
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
               else if(hasReceivedStatusFromParent[i]==true || nbr[i].getDepth()==min_depth){
                   nbr_exists[i]=true;
                   nbr_isRefined[i]=false;
               }
               else{// the neighbor is going to get derefined
                   nbr_exists[i]=false;
               }
	   }
	   else if(nbr_exists[i] && nbr_isRefined[i] && nbr_decision[i]!=STAY){
	      nbr_isRefined[i]=false;
	   }
       }
    }
    else if(decision==REFINE){//I will Now become Inactive and therefore I need not Store Neighbor Status
    }
    if(decision==STAY || (isRefined && !isGrandParent && !parentHasAlreadyMadeDecision))
        CkStartQD(CkIndex_Advection::doStep(), &thishandle);
}

void Advection::recvChildData(ChildDataMsg *msg){
    if(!hasAllocatedMemory){
        hasAllocatedMemory=true;
	mem_allocate_all();
    }
    int st_i, end_i, st_j, end_j;
    if(msg->childNum==0){
    	st_i=block_width/2+1; end_i=block_width;
	st_j=1; end_j=block_height;
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
        ckout << "Error: recvChildData(ChildDataMsg*) received " << msg->childNum << "as ChildNum" <<endl;
	CkExit();
    }

    int ctr=0;
    for(int i=st_i; i<=end_i; i++){
        for(int j=st_j; j<=end_j; j++){
	    u[index(i,j)]=msg->child_u[ctr++];
	}
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
    
    setNbrStatus(c1, msg);
    setNbrStatus(c2, msg);

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
        else if(hasReceivedStatusFromParent[dir]==true){
            nbr_exists[dir]=true;
            nbr_isRefined[dir]=false;
        }
        else{
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
           // ckout << u[i+1][j+1] << endl;
            intensity[3*(i+j*w)+ 0] = 255;// red component
            intensity[3*(i+j*w)+ 1] = 255 - (u[index(i+1,j+1)]-1)*255;// BLUE component
            intensity[3*(i+j*w)+ 2] = 255 - (u[index(i+1,j+1)]-1)*255;// green component
        }
    }
    liveVizDeposit(m, sy*w,sx*h, w,h, intensity, this);
    delete[] intensity;
}
#define index_l(i,j)  (int)((j)*(block_width) + i)

void Advection::interpolate(double *u, double *refined_u, int xstart, int xend, int ystart, int yend){
    double sx, sy;
    int m=1, n=1;
    for(int i=xstart; i<=xend; i++){
        for(int j=ystart; j<=yend; j++){
            sx = (u[index(i-1,j)]-u[index(i+1,j)])/(2*dx);
            sy = (u[index(i,j-1)]-u[index(i,j+1)])/(2*dy);

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
    double sx, sy, s;
    size_t sz = (block_width)*(block_height); 
    double *refined_u = new double[sz];

    interpolate(u, refined_u, 1, block_width/2, 1, block_height/2);
    //initialize the child
    InitRefineMsg * msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
    thisProxy(thisIndex.getChild("01")).insert(msg);

    interpolate(u, refined_u, block_width/2+1, block_width, 1, block_height/2);
    //initialize the child
    msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
    thisProxy(thisIndex.getChild("00")).insert(msg);

    interpolate(u, refined_u, 1, block_width/2, block_height/2+1, block_height);
    //init the child
    msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
    thisProxy(thisIndex.getChild("10")).insert(msg);

    interpolate(u, refined_u, block_width/2+1, block_width, block_height/2+1, block_height);
    //init the child
    msg = new (sz, NUM_NEIGHBORS, NUM_NEIGHBORS, 3*NUM_NEIGHBORS)InitRefineMsg(dx/2, dy/2, myt, mydt, iterations, refined_u, nbr_exists, nbr_isRefined, nbr_decision);
    thisProxy(thisIndex.getChild("11")).insert(msg);
    thisProxy.doneInserting();
    //delete [] refined_u;
    ckout << thisIndex.getIndexString() << " done with refinement" << endl;
}

Advection::Advection(InitRefineMsg* msg){
    __sdag_init();
    usesAtSync=CmiTrue;

//Called as a result of refinement of parent
    ckout << "Inserting New Zone: " << thisIndex.getIndexString() << endl;
    CBase_Advection();
    this->exists = true;
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
                    if(msg->parent_nbr_decision[getNbrDir(thisIndex.getQuadI(), dir)]==REFINE)
                        nbr_isRefined[dir]=true;
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
    hasReceived = *new set<int>();
    dx = msg->dx;
    dy = msg->dy;

    myt = msg->myt;
    mydt = msg->mydt;
    iterations = msg->iterations;

    mem_allocate_all();

    //x[i] and y[i] need not be initialized they are needed only for setting initial conditions
    delete [] x;
    delete [] y;
 	
    //Initialize u - For boundaries I have to wait for the neighbors
    //to send the values, rest of it can be initialized by the values 
    //received from the parent
    int ctr=0;
    for(int j=1; j<=block_height; j++)
    	for(int i=1; i<=block_width; i++)
	    u[index(i,j)]=msg->refined_u[ctr++];

    
    //delete the message
    delete msg;

    //call Quiesence detection to begin next iteration*/
    CkStartQD(CkIndex_Advection::doStep(), &thishandle);
}
#include "Advection.def.h"
