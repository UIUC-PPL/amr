#ifndef __CONSTANTS__
#define __CONSTANTS__

#define NUM_NEIGHBORS 4
#define NUM_CHILDREN 4

#define FOR_EACH_ZONE for ( int i=1; i<=block_width; i++) for(int j=1; j<=block_height; j++){
#define END_FOR }
#define FOR_EACH_NEIGHBOR for( int i=0; i<NUM_NEIGHBORS; i++){
#define FOR_EACH_CHILD for( int i=0; i<NUM_CHILDREN; i++){
#endif
