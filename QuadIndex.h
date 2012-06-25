#include <string>

const int bits_per_int = 8*sizeof(int);

inline bool liesin(int x, int st, int end){
    if (x<=end && x>=st)
        return true;
    else return false;
}

enum DIR {  UP=0, DOWN=1, LEFT=2, RIGHT=3, LEFT_UP=4, LEFT_DOWN=5, RIGHT_UP=6, 
            RIGHT_DOWN=7, UP_LEFT=8, UP_RIGHT=9, DOWN_LEFT=10, DOWN_RIGHT=11};
PUPbytes(DIR);

static DIR SENDER_DIR[NUM_NEIGHBORS] = {DOWN,UP,RIGHT,LEFT};

enum DECISION {INV=-1, DEREFINE=0, STAY=1, REFINE=2};
PUPbytes(DECISION);


class QuadIndex{
public:
    int nbits;
    unsigned int bitVector;// bit vector storing the index bits
    // number of bits in the index
        
public:
    QuadIndex(){bitVector=0;nbits=0;}
    QuadIndex(const QuadIndex& qindex){bitVector = qindex.bitVector; nbits = qindex.nbits;}
    QuadIndex(int bitVector, int nbits){
        this->bitVector = bitVector; 
        this->nbits = nbits;
    }
    QuadIndex(const char*);
    QuadIndex(int, int, int);

    void operator = (const QuadIndex& qidx) {
        nbits = qidx.nbits;
        bitVector = qidx.bitVector;
    }

    bool operator !=(QuadIndex qidx) const{
        if(nbits==qidx.nbits && bitVector==qidx.bitVector)
            return false;
        return true;
    }
    
    bool operator ==(QuadIndex qidx) const{
        if(nbits==qidx.nbits && getIndexString() == qidx.getIndexString())
            return true;
        return false;
    }

    void getCoordinates(int&, int&)const;
    void getPhysicalCoordinates(int &x, int &y) const;
    std::string getIndexString() const;
    inline int getDepth() const{
        return nbits/2;
    }

    QuadIndex getNeighbor(int) const;
    QuadIndex getParent() const;

    QuadIndex getChild(unsigned int) const;
    int getChildNum() const;
    int getQuadI() const;
    
    unsigned int getBitVector() const{ return bitVector;}
    unsigned int getNbits() const{ return nbits;}

    void pup(PUP::er &p){
        p|nbits;
        p|bitVector;
    }

    void getSiblingInDirection(DIR, int&, int&) const;

    bool operator<(const QuadIndex& rhs) const {
      if (nbits < rhs.nbits) return true;
      if (nbits > rhs.nbits) return false;
      if (bitVector < rhs.bitVector) return true;
      return false;
    }
};

class CkArrayIndexQuadIndex: public CkArrayIndex {
public:
    QuadIndex *idx;
    
public:
    CkArrayIndexQuadIndex(const QuadIndex &in){
        idx = new(index) QuadIndex(in);
        nInts = sizeof(in)/sizeof(int);
    }
    CkArrayIndexQuadIndex(){
        idx = new(index) QuadIndex();
        nInts = sizeof(QuadIndex)/sizeof(int);
    }
    
    CkArrayIndexQuadIndex &operator=(const CkArrayIndexQuadIndex &that)  {
        nInts = that.nInts;
        memcpy(idx, that.idx, sizeof(int)*2);
        return *this;
    }
};
