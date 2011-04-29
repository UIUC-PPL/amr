#include <cmath>
#include <cstring>
#include "math.h"
#include "charm++.h"
#include "Constants.h"
#include <map>
#include "boost/assign.hpp"
#include "QuadIndex.h"

using namespace std;
using namespace boost::assign;

extern map<char*, DIR> nbrDirectionMap;
extern map<DIR, DIR> reverse_dir_map;

/*constructor
__________
|    |    |
|_01_|_00_|        
|    |    |
|_10_|_11_|

This is how quadrants are defined
*/
QuadIndex::QuadIndex(const char* index){
    nbits = strlen(index);
    bitVector=0;
    //cout << "nbits: " << nbits << endl;
    int nElements = ceil((double)nbits/bits_per_int);// number of integers requries in the bit vector
    
    if(nElements > 1){
        ckout << "Tree Size exceeded" << endl;
        return;
    }
    if(nbits%2!=0){
        ckout << "Invalid Index String, Index length should be a multiple of 2" << endl;
        ckout << "nbits: " << nbits << ", " << index << endl;
        return;
    }
    //reset bit vector
    
    bitVector = 0;
    
    //set bit vector
    for(int i=0; i<nbits; i++){
        int bit = (int)index[i]-48;
        
        if(bit!=0 && bit !=1){
            ckout << "Invalid Index String" <<endl;
            return;
        }
        bitVector = bitVector | bit<<(bits_per_int - 1 - i);
    }
    //cout << bitVector << endl;
}
    
    //contructor- given the coordinates of the node and its depth, generate the index vector
QuadIndex::QuadIndex(int x, int y, int depth){ 
    bitVector = 0; 
    nbits = 2*depth;
    int quad;
    int range = std::pow(2.0, (double)depth);
    int r1 = 0, r2 = range-1, c1 = 0, c2 = range-1;
    
    int index = 0;
    for(int i=0; i<depth; i++){
        if(liesin(y, 0, (r1+r2)/2) && liesin(x, (c1+c2)/2+1, c2))quad=0;
        else if (liesin(y, 0, (r1+r2)/2) && liesin(x, c1, (c1+c2)/2))quad=1;
        else if (liesin(y, (r1+r2)/2 + 1, r2) && liesin(x, c1, (c1+c2)/2))quad=2;
        else quad=3;
        
        //set bitVector
        int bit = quad/2;
        bitVector = bitVector | bit<<(bits_per_int - 1 - index);
        bit = quad%2;
        bitVector = bitVector | bit<<(bits_per_int - 1 - (index+1));
        index += 2;
         
         r2 = r2 - std::pow(2.0, (double)depth-i-1);
         c2 = c2 - std::pow(2.0, (double)depth-i-1);
         
         x =  x%(r2+1);
         y =  y%(r2+1);
    }
}

void QuadIndex::getCoordinates(int &x, int &y) const{
    int depth = nbits/2;
    int r1 = 0, r2 = std::pow(2.0, (double)depth)-1, c1 = 0, c2 = std::pow(2.0, (double)depth)-1;
    for(int i=0; i<nbits; i+=2){
        int bit0 = (bitVector & (1<<(bits_per_int - 1 - i)))>0?1:0;
        int bit1 = (bitVector & (1<<(bits_per_int - 1 - (i+1))))>0?1:0;
        
        int quad = 2*bit0 + bit1;
     //   cout << "quad: " << quad << endl;
        if(quad==0){
            r2 = (r1+r2)/2;
            c1 = (c1+c2)/2+1;
        }
        else if(quad==1){
            r2 = (r1+r2)/2;
            c2 = (c1+c2)/2;
        }
        else if(quad==2){
            r1 = (r1+r2)/2 + 1;
            c2 = (c1+c2)/2;
        }
        else{
            r1 = (r1+r2)/2 + 1;
            c1 = (c1+c2)/2 + 1;
        }
//           cout << r1 << ", " << r2 << ", " << c1 << ", " << c2 << endl;
    }
//       cout << "getCoordinates for " << getIndexString() << ": " << x << ", " << y << endl;     
    x = c1;//at the end of for loop r1=r2
    y = r1;// at the end of for loop c1=c2
}
    
    // returns the Index String
char* QuadIndex::getIndexString() const{
    char* str = new char[nbits+1];
    int i;
    for(i=0; i<nbits; i++){
        str[i] = (bitVector & (1<<(bits_per_int - 1 - i)))>0?49:48;
    }
    str[i]='\0';
    return str;
}
    
    
QuadIndex QuadIndex::getNeighbor(int dir) const{
    if(dir==UP){
        int x, y;
        getCoordinates(x, y);
        int depth = nbits/2;
        int range = std::pow(2.0, (double)depth);
        int yc = (y==0)?(range-1):y-1;
        int xc = x;
        return *new QuadIndex(xc, yc, depth);
    }
    else if(dir==DOWN){
        int x, y;
        getCoordinates(x, y);
        int depth = nbits/2;
        int range = std::pow(2.0, (double)depth);
        int yc = (y+1)%range;
        int xc = x;
        return *new QuadIndex(xc, yc, depth);
    }
    else if(dir==LEFT){
        int x, y;
        getCoordinates(x, y);
        int depth = nbits/2;
        int range = std::pow(2.0, (double)depth);
        int yc = y;
        int xc = (x==0)?(range-1):x-1;
        return *new QuadIndex(xc, yc, depth);
    }
    else if(dir==RIGHT){
       int x, y;
        getCoordinates(x, y);
        int depth = nbits/2;
        int range = std::pow(2.0, (double)depth);
        int yc = y;
        int xc = (x+1)%range;
        return *new QuadIndex(xc, yc, depth);
    }
}
    
QuadIndex QuadIndex::getParent() const{
    char* myIndexString = new char[nbits+1];
    memcpy(myIndexString, getIndexString(), nbits);
    myIndexString[nbits-2]='\0';//Make it a parent Index String
    char* parentIndex = myIndexString;

    return *new QuadIndex(parentIndex);
}

QuadIndex QuadIndex::getChild(const char* child_num) const{
    char* str = new char[nbits+3];
    str[nbits+2]='\0';
    memcpy(str, getIndexString(), nbits);
    str[nbits]='\0';
    strcat(str, child_num);
    return *new QuadIndex(str);
}
    
QuadIndex QuadIndex::getChild(int idx) const{
    switch(idx){
        case 0:
            return getChild("00");
            break;
        case 1:
            return getChild("01");
            break;
        case 2:
            return getChild("10");
            break;
        case 3:
            return getChild("11");
            break;
        default:
            CkAbort("Error in getChild()");
            break;
    }
}

int QuadIndex::getChildNum() const{
    if (nbits==0){//When I am the root
        ckout << "CAUTION: getChildNum() called on Root of tree" << endl;
        return -1;
    }
    int bit0 = ((bitVector & 1<<(bits_per_int-nbits))>0)?1:0;
    int bit1 = ((bitVector & 1<<(bits_per_int-nbits+1))>0)?1:0;
    return 2*bit1+bit0;
}

char* QuadIndex::getQuadC() const{
    char* index = getIndexString();
    char* quad = new char[3];
    int len = strlen(index);
    memcpy(quad, index + len -2, 2);
    quad[2] = 0;
    delete [] index;
    return quad;
}

int QuadIndex::getQuadI() const{
    char* index = getIndexString();
    char* quad = new char[3];
    int len = strlen(index);
    memcpy(quad, index + len -2, 2);
    quad[2] = 0;
    delete [] index;
    if(strcmp(quad, "00")==0) return 0;
    else if(strcmp(quad, "01")==0)return 1;
    else if(strcmp(quad, "10")==0)return 2;
    else return 3;
}


DIR QuadIndex::getSiblingDirection(QuadIndex nbr) const{
    char* quad1 = getQuadC();
    char* quad2 = nbr.getQuadC();

    char* key = strcat(quad1, quad2);
    delete [] quad2;
    DIR dir;
    if(this->getParent() == nbr.getParent())
        dir = nbrDirectionMap.find(key)->second;
    else
        dir= reverse_dir_map[nbrDirectionMap.find(key)->second];

    delete [] key;
    return dir;
}

void QuadIndex::getSiblingInDirection(DIR dir, int &c1, int &c2) const{
    if(dir==RIGHT){
        c1=0; c2=3;}
    else if(dir==UP){
        c1=0; c2=1;}
    else if(dir==LEFT){
        c1=1; c2=2;}
    else if(dir==DOWN){
        c1=2; c2=3;}
}


