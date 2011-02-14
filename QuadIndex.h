
const int bits_per_int = 8*sizeof(int);

inline bool liesin(int x, int st, int end){
    if (x<=end && x>=st)
        return true;
    else return false;
}

enum DIR {  UP=0, DOWN=1, LEFT=2, RIGHT=3, LEFT_UP, LEFT_DOWN, RIGHT_UP, 
            RIGHT_DOWN, UP_LEFT, UP_RIGHT, DOWN_LEFT, DOWN_RIGHT};

static int SENDER_DIR[NUM_NEIGHBORS] = {1,0,3,2};

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

    void getCoordinates(int&, int&)const;
    char* getIndexString() const;

    QuadIndex getNeighbor(int) const;
    QuadIndex getParent() const;
    QuadIndex getChild(const char*)const;

    QuadIndex getChild(int) const;
    int getChildNum() const;
    
    unsigned int getBitVector() const{ return bitVector;}
    unsigned int getNbits() const{ return nbits;}

    void pup(PUP::er &p){
        p|nbits;
        p|bitVector;
    }
};

class CkArrayIndexQuadIndex: public CkArrayIndex {
private:
    QuadIndex *idx;
    
public:
    CkArrayIndexQuadIndex(const QuadIndex &in){
        idx = new(index) QuadIndex(in);
        nInts = sizeof(in)/sizeof(int);
    }
    CkArrayIndexQuadIndex(){
        idx = new(index) QuadIndex();
    }
    
    CkArrayIndexQuadIndex &operator=(const CkArrayIndexQuadIndex &that)  {
        nInts = that.nInts;
        memcpy(idx, that.idx, sizeof(int)*2);
        return *this;
    }
};
