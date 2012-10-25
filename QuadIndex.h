#ifndef QUADINDEX_H
#define QUADINDEX_H

using std::string;

class QuadIndex{
public:
  int nbits;// number of bits in the index
  unsigned int bitVector;// bit vector storing the index bits
        
public:
  QuadIndex():bitVector(0), nbits(0){}
  QuadIndex(const QuadIndex& qindex):bitVector(qindex.bitVector), nbits(qindex.nbits){}
  QuadIndex(int _bitVector, int _nbits):bitVector(_bitVector), nbits(_nbits){}
  QuadIndex(int x, int y, int depth);

  void operator = (const QuadIndex& qidx) {
    nbits = qidx.nbits;
    bitVector = qidx.bitVector;
  }

  bool operator != (QuadIndex qidx) const { return !(*this == qidx); }
    
  bool operator == (QuadIndex qidx) const {
    return (nbits==qidx.nbits && getIndexString() == qidx.getIndexString())? true: false;
  }

  void getCoordinates(int&, int&)const;
  std::string getIndexString() const;
  inline int getDepth() const{ return nbits/2; }

  QuadIndex getNeighbor(int) const;
  QuadIndex getParent() const;

  QuadIndex getChild(unsigned int) const;
  int getQuadrant() const;
    
  unsigned int getBitVector() const{ return bitVector;}
  unsigned int getNbits() const{ return nbits;}

  void pup(PUP::er &p){
    p|nbits;
    p|bitVector;
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

#endif
