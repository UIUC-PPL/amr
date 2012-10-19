#ifndef __CONSTANTS__
#define __CONSTANTS__

#ifndef AMR_REVISION
#define AMR_REVISION Unknown
#endif
#define _QUOTEIT(x) #x

#define INQUOTES(x) _QUOTEIT(x)
#ifdef LOGGER
#define VB(x) do { x } while(false)
#else
#define VB(x) do { } while(false)
#endif

#define NUM_NEIGHBORS 4
#define NUM_CHILDREN 4

#define FOR_EACH_ZONE for ( int i=1; i<=block_width; i++) for(int j=1; j<=block_height; j++){
#define FOR_EACH_NEIGHBOR for( int i=0; i<NUM_NEIGHBORS; i++){
#define FOR_EACH_CHILD for( int i=0; i<NUM_CHILDREN; i++){
#define END_FOR }

#define index(i,j)      (int)((j)*(block_width+2) + i)
#define index_l(i,j)    (int)((j)*(block_width) + i)
#define index_c(i,j)    (int)((j)*(block_width/2) + i)
#define liesin(x, st, end) ((x<=end && x>=st) ? true:false)

#define isLeaf !isRefined

const int bits_per_int = 8*sizeof(int);

enum Dir {  UP=0, DOWN=1, LEFT=2, RIGHT=3, LEFT_UP=4, LEFT_DOWN=5, RIGHT_UP=6, 
            RIGHT_DOWN=7, UP_LEFT=8, UP_RIGHT=9, DOWN_LEFT=10, DOWN_RIGHT=11};
PUPbytes(Dir);

enum Decision {INV=-1, COARSEN=0, STAY=1, REFINE=2};
PUPbytes(Decision);

#endif
