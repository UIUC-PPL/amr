#ifndef _DECL_Advection_H_
#define _DECL_Advection_H_
#include "charm++.h"
/* DECLS: array Advection: ArrayElement{
Advection(CkMigrateMessage* impl_msg);
Advection(void);
Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
void refine(void);
void inform_nbr_of_refinement(void);
void receiveGhosts(int iter, int dir, int width, const double *u);
void begin_iteration(void);
void doStep(void);
void requestNextFrame(liveVizRequestMsg* impl_msg);
void setReal(void);
void setDummy(void);
};
 */
 class Advection;
 class CkIndex_Advection;
 class CProxy_Advection;
 class CProxyElement_Advection;
 class CProxySection_Advection;
/* --------------- index object ------------------ */
class CkIndex_Advection:public CProxyElement_ArrayElement{
  public:
    typedef Advection local_t;
    typedef CkIndex_Advection index_t;
    typedef CProxy_Advection proxy_t;
    typedef CProxyElement_Advection element_t;
    typedef CProxySection_Advection section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Advection(CkMigrateMessage* impl_msg);
 */
    static int __idx_Advection_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_Advection_CkMigrateMessage; }
    static void _call_Advection_CkMigrateMessage(void* impl_msg,Advection* impl_obj);

/* DECLS: Advection(void);
 */
    static int __idx_Advection_void;
    static int ckNew(void) { return __idx_Advection_void; }
    static void _call_Advection_void(void* impl_msg,Advection* impl_obj);

/* DECLS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */
    static int __idx_Advection_marshall2;
    static int ckNew(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2) { return __idx_Advection_marshall2; }
    static void _call_Advection_marshall2(void* impl_msg,Advection* impl_obj);
    static int _callmarshall_Advection_marshall2(char* impl_buf,Advection* impl_obj);
    static void _marshallmessagepup_Advection_marshall2(PUP::er &p,void *msg);

/* DECLS: void refine(void);
 */
    static int __idx_refine_void;
    static int refine(void) { return __idx_refine_void; }
    static void _call_refine_void(void* impl_msg,Advection* impl_obj);

/* DECLS: void inform_nbr_of_refinement(void);
 */
    static int __idx_inform_nbr_of_refinement_void;
    static int inform_nbr_of_refinement(void) { return __idx_inform_nbr_of_refinement_void; }
    static void _call_inform_nbr_of_refinement_void(void* impl_msg,Advection* impl_obj);

/* DECLS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
    static int __idx_receiveGhosts_marshall5;
    static int receiveGhosts(int iter, int dir, int width, const double *u) { return __idx_receiveGhosts_marshall5; }
    static void _call_receiveGhosts_marshall5(void* impl_msg,Advection* impl_obj);
    static int _callmarshall_receiveGhosts_marshall5(char* impl_buf,Advection* impl_obj);
    static void _marshallmessagepup_receiveGhosts_marshall5(PUP::er &p,void *msg);

/* DECLS: void begin_iteration(void);
 */
    static int __idx_begin_iteration_void;
    static int begin_iteration(void) { return __idx_begin_iteration_void; }
    static void _call_begin_iteration_void(void* impl_msg,Advection* impl_obj);

/* DECLS: void doStep(void);
 */
    static int __idx_doStep_void;
    static int doStep(void) { return __idx_doStep_void; }
    static void _call_doStep_void(void* impl_msg,Advection* impl_obj);

/* DECLS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
    static int __idx_requestNextFrame_liveVizRequestMsg;
    static int requestNextFrame(liveVizRequestMsg* impl_msg) { return __idx_requestNextFrame_liveVizRequestMsg; }
    static void _call_requestNextFrame_liveVizRequestMsg(void* impl_msg,Advection* impl_obj);

/* DECLS: void setReal(void);
 */
    static int __idx_setReal_void;
    static int setReal(void) { return __idx_setReal_void; }
    static void _call_setReal_void(void* impl_msg,Advection* impl_obj);

/* DECLS: void setDummy(void);
 */
    static int __idx_setDummy_void;
    static int setDummy(void) { return __idx_setDummy_void; }
    static void _call_setDummy_void(void* impl_msg,Advection* impl_obj);

};
/* --------------- element proxy ------------------ */
 class CProxyElement_Advection : public CProxyElement_ArrayElement{
  public:
    typedef Advection local_t;
    typedef CkIndex_Advection index_t;
    typedef CProxy_Advection proxy_t;
    typedef CProxyElement_Advection element_t;
    typedef CProxySection_Advection section_t;

    CProxyElement_Advection(void) {}
    CProxyElement_Advection(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxyElement_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxyElement_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxyElement_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxyElement_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxyElement_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxyElement_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxyElement_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxyElement_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxyElement_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxyElement_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxyElement_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxyElement_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxyElement_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxyElement_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxyElement_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxyElement_ArrayElement::ckSetReductionClient(cb); }
inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
  { CProxyElement_ArrayElement::ckInsert(m,ctor,onPe); }
inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
  { CProxyElement_ArrayElement::ckSend(m,ep,opts); }
inline void *ckSendSync(CkArrayMessage *m, int ep) const
  { return CProxyElement_ArrayElement::ckSendSync(m,ep); }
inline const CkArrayIndex &ckGetIndex() const
  { return CProxyElement_ArrayElement::ckGetIndex(); }
    Advection *ckLocal(void) const
      { return (Advection *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_Advection(const CkArrayID &aid,const CkArrayIndexQuadIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS) {}
    CProxyElement_Advection(const CkArrayID &aid,const CkArrayIndexQuadIndex &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: Advection(CkMigrateMessage* impl_msg);
 */

/* DECLS: Advection(void);
 */
    void insert(int onPE=-1);
/* DECLS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */
    void insert(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void refine(void);
 */
    void refine(void) ;

/* DECLS: void inform_nbr_of_refinement(void);
 */
    void inform_nbr_of_refinement(void) ;

/* DECLS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
    void receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void begin_iteration(void);
 */
    void begin_iteration(void) ;

/* DECLS: void doStep(void);
 */
    void doStep(void) ;

/* DECLS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
    void requestNextFrame(liveVizRequestMsg* impl_msg) ;

/* DECLS: void setReal(void);
 */
    void setReal(void) ;

/* DECLS: void setDummy(void);
 */
    void setDummy(void) ;

};
PUPmarshall(CProxyElement_Advection)
/* ---------------- collective proxy -------------- */
 class CProxy_Advection : public CProxy_ArrayElement{
  public:
    typedef Advection local_t;
    typedef CkIndex_Advection index_t;
    typedef CProxy_Advection proxy_t;
    typedef CProxyElement_Advection element_t;
    typedef CProxySection_Advection section_t;

    CProxy_Advection(void) {}
    CProxy_Advection(const ArrayElement *e) : CProxy_ArrayElement(e){  }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxy_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxy_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxy_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxy_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxy_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxy_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxy_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxy_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxy_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxy_ArrayElement::ckSetReductionClient(cb); }
    static CkArrayID ckNew(void) {return ckCreateEmptyArray();}
//Generalized array indexing:
    CProxyElement_Advection operator [] (const CkArrayIndexQuadIndex &idx) const
        {return CProxyElement_Advection(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_Advection operator() (const CkArrayIndexQuadIndex &idx) const
        {return CProxyElement_Advection(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxy_Advection(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_Advection(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: Advection(CkMigrateMessage* impl_msg);
 */

/* DECLS: Advection(void);
 */
    static CkArrayID ckNew(const CkArrayOptions &opts);

/* DECLS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */
    static CkArrayID ckNew(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void refine(void);
 */
    void refine(void) ;

/* DECLS: void inform_nbr_of_refinement(void);
 */
    void inform_nbr_of_refinement(void) ;

/* DECLS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
    void receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void begin_iteration(void);
 */
    void begin_iteration(void) ;

/* DECLS: void doStep(void);
 */
    void doStep(void) ;

/* DECLS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
    void requestNextFrame(liveVizRequestMsg* impl_msg) ;

/* DECLS: void setReal(void);
 */
    void setReal(void) ;

/* DECLS: void setDummy(void);
 */
    void setDummy(void) ;

};
PUPmarshall(CProxy_Advection)
/* ---------------- section proxy -------------- */
 class CProxySection_Advection : public CProxySection_ArrayElement{
  public:
    typedef Advection local_t;
    typedef CkIndex_Advection index_t;
    typedef CProxy_Advection proxy_t;
    typedef CProxyElement_Advection element_t;
    typedef CProxySection_Advection section_t;

    CProxySection_Advection(void) {}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxySection_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxySection_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxySection_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxySection_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxySection_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxySection_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxySection_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxySection_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxySection_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxySection_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxySection_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxySection_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxySection_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxySection_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxySection_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxySection_ArrayElement::ckSetReductionClient(cb); }
inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
 { CProxySection_ArrayElement::ckSend(m,ep,opts); }
inline CkSectionInfo &ckGetSectionInfo()
  { return CProxySection_ArrayElement::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
  { return CProxySection_ArrayElement::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
  { return CProxySection_ArrayElement::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
  { return CProxySection_ArrayElement::ckGetSectionID(i); }
inline CkArrayID ckGetArrayIDn(int i) const
{return CProxySection_ArrayElement::ckGetArrayIDn(i); } 
inline CkArrayIndexMax *ckGetArrayElements() const
  { return CProxySection_ArrayElement::ckGetArrayElements(); }
inline CkArrayIndexMax *ckGetArrayElements(int i) const
{return CProxySection_ArrayElement::ckGetArrayElements(i); }
inline int ckGetNumElements() const
  { return CProxySection_ArrayElement::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{return CProxySection_ArrayElement::ckGetNumElements(i); } 
//Generalized array indexing:
    CProxyElement_Advection operator [] (const CkArrayIndexQuadIndex &idx) const
        {return CProxyElement_Advection(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_Advection operator() (const CkArrayIndexQuadIndex &idx) const
        {return CProxyElement_Advection(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_Advection(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_Advection(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_Advection(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    CProxySection_Advection(int n, const CkArrayID *aid, CkArrayIndexMax const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_Advection(int n, const CkArrayID *aid, CkArrayIndexMax const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: Advection(CkMigrateMessage* impl_msg);
 */

/* DECLS: Advection(void);
 */

/* DECLS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */

/* DECLS: void refine(void);
 */
    void refine(void) ;

/* DECLS: void inform_nbr_of_refinement(void);
 */
    void inform_nbr_of_refinement(void) ;

/* DECLS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
    void receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void begin_iteration(void);
 */
    void begin_iteration(void) ;

/* DECLS: void doStep(void);
 */
    void doStep(void) ;

/* DECLS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
    void requestNextFrame(liveVizRequestMsg* impl_msg) ;

/* DECLS: void setReal(void);
 */
    void setReal(void) ;

/* DECLS: void setDummy(void);
 */
    void setDummy(void) ;

};
PUPmarshall(CProxySection_Advection)
#define Advection_SDAG_CODE                                                    \
public:                                                                        \
  void doStep() {                                                              \
    _slist_0();                                                                \
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _dummyEP, CkMyPe(), 0, NULL);        \
  }                                                                            \
                                                                               \
private:                                                                       \
  void doStep_end() {                                                          \
  }                                                                            \
                                                                               \
  void _slist_0() {                                                            \
    _atomic_0();                                                               \
  }                                                                            \
  void _slist_0_end() {                                                        \
    doStep_end();                                                              \
  }                                                                            \
  void _atomic_0() {                                                           \
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, __idx_Advection_begin_iteration, CkMyPe(), 0, NULL); \
                                                                               \
                    begin_iteration();                                         \
                                                                               \
    _TRACE_END_EXECUTE();                                                      \
    _for_0();                                                                  \
  }                                                                            \
                                                                               \
  void _for_0() {                                                              \
    imsg = 0;                                                                  \
    if ( imsg < 4) {                                                           \
      _slist_1();                                                              \
    } else {                                                                   \
      _atomic_2();                                                             \
    }                                                                          \
  }                                                                            \
  void _for_0_end() {                                                          \
 imsg++;                                                                       \
    if ( imsg < 4) {                                                           \
      _slist_1();                                                              \
    } else {                                                                   \
      _atomic_2();                                                             \
    }                                                                          \
  }                                                                            \
  void _slist_1() {                                                            \
    _when_0();                                                                 \
  }                                                                            \
  void _slist_1_end() {                                                        \
    _for_0_end();                                                              \
  }                                                                            \
  int _when_0() {                                                              \
    CMsgBuffer *receiveGhosts_buf;                                             \
    CkMarshallMsg *receiveGhosts_msg;                                          \
                                                                               \
    receiveGhosts_buf = __cDep->getMessage(0, iterations);                     \
                                                                               \
    if ((receiveGhosts_buf != 0)) {                                            \
       receiveGhosts_msg = (CkMarshallMsg *)receiveGhosts_buf->msg;            \
       char *receiveGhosts_impl_buf=((CkMarshallMsg *)receiveGhosts_msg)->msgBuf;\
       PUP::fromMem receiveGhosts_implP(receiveGhosts_impl_buf);               \
       int iter; receiveGhosts_implP|iter;                                     \
       int dir; receiveGhosts_implP|dir;                                       \
       int width; receiveGhosts_implP|width;                                   \
      int impl_off_u; receiveGhosts_implP|impl_off_u;                          \
       receiveGhosts_impl_buf+=CK_ALIGN(receiveGhosts_implP.size(),16);        \
    double *u=(double *)(receiveGhosts_impl_buf+impl_off_u);                   \
       __cDep->removeMessage(receiveGhosts_buf);                               \
       delete receiveGhosts_buf;                                               \
       _atomic_1(iter, dir, width, u);                                         \
       delete receiveGhosts_msg;                                               \
       return 1;                                                               \
    } else {                                                                   \
       CWhenTrigger *tr;                                                       \
       tr = new CWhenTrigger(0, 1, 1, 0);                                      \
       tr->args[0] = (size_t) CkAllocSysMsg();                                 \
       tr->entries[0] = 0;                                                     \
       tr->refnums[0] = iterations;                                            \
       __cDep->Register(tr);                                                   \
       return 0;                                                               \
    }                                                                          \
  }                                                                            \
                                                                               \
  void _when_0_end(int iter, int dir, int width, double *u) {                  \
    _slist_1_end();                                                            \
  }                                                                            \
                                                                               \
  void _atomic_1(int iter, int dir, int width, double *u) {                    \
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, __idx_Advection_process_ghosts, CkMyPe(), 0, NULL); \
                                                                               \
                        process(dir, width, u);                                \
                                                                               \
    _TRACE_END_EXECUTE();                                                      \
    _when_0_end(iter, dir, width, u);                                          \
  }                                                                            \
                                                                               \
  void _atomic_2() {                                                           \
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, __idx_Advection_doWork, CkMyPe(), 0, NULL); \
                                                                               \
                    compute_and_iterate();                                     \
                                                                               \
    _TRACE_END_EXECUTE();                                                      \
    _slist_0_end();                                                            \
  }                                                                            \
                                                                               \
public:                                                                        \
  void receiveGhosts(int iter, int dir, int width, double *u) {                \
    CWhenTrigger *tr;                                                          \
    void* _bgParentLog = NULL;                                                 \
    CMsgBuffer* cmsgbuf;                                                       \
    int impl_off=0; int impl_arrstart=0;                                       \
    int impl_off_u, impl_cnt_u;                                                \
    impl_off_u=impl_off=CK_ALIGN(impl_off,sizeof(double));                     \
    impl_off+=(impl_cnt_u=sizeof(double)*(width));                             \
    {                                                                          \
      PUP::sizer implP1;                                                       \
      implP1|iter;                                                             \
      implP1|dir;                                                              \
      implP1|width;                                                            \
      implP1|impl_off_u;                                                       \
      impl_arrstart=CK_ALIGN(implP1.size(),16);                                \
      impl_off+=impl_arrstart;                                                 \
    }                                                                          \
    CkMarshallMsg *impl_msg1=CkAllocateMarshallMsg(impl_off,NULL);             \
    {                                                                          \
      PUP::toMem implP1((void *)impl_msg1->msgBuf);                            \
      implP1|iter;                                                             \
      implP1|dir;                                                              \
      implP1|width;                                                            \
      implP1|impl_off_u;                                                       \
    }                                                                          \
    char *impl_buf1=impl_msg1->msgBuf+impl_arrstart;                           \
    memcpy(impl_buf1+impl_off_u,u,impl_cnt_u);                                 \
   cmsgbuf = __cDep->bufferMessage(0,(void *) impl_msg1, (void*) _bgParentLog,iter);\
    tr = __cDep->getTrigger(0,iter);                                           \
    if (tr == 0)                                                               \
      return;                                                                  \
    _TRACE_END_EXECUTE();                                                      \
    CkFreeSysMsg((void  *)tr->args[0]);                                        \
    _when_0();                                                                 \
    delete tr;                                                                 \
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _dummyEP, CkMyPe(), 0, NULL);        \
    return;                                                                    \
  }                                                                            \
                                                                               \
private:                                                                       \
  CDep *__cDep;                                                                \
  void __sdag_init(void) {                                                     \
    __cDep = new CDep(1,1);                                                    \
    __cDep->addDepends(0,0);                                                   \
  }                                                                            \
public:                                                                        \
  void __sdag_pup(PUP::er& p) {                                                \
    if (__cDep) { __cDep->pup(p); }                                            \
  }                                                                            \
  static void __sdag_register() {                                              \
                                                                               \
    __idx_Advection_begin_iteration = CkRegisterEp("Advection_begin_iteration(void)", NULL, 0, CkIndex_Advection::__idx, 0);\
    __idx_Advection_process_ghosts = CkRegisterEp("Advection_process_ghosts(void)", NULL, 0, CkIndex_Advection::__idx, 0);\
    __idx_Advection_doWork = CkRegisterEp("Advection_doWork(void)", NULL, 0, CkIndex_Advection::__idx, 0);\
  }                                                                            \
  static int __idx_Advection_begin_iteration;                                  \
  static int __idx_Advection_process_ghosts;                                   \
  static int __idx_Advection_doWork;                                           \

#define Advection_SDAG_CODE_DEF \
  int Advection::__idx_Advection_begin_iteration=0;\
  int Advection::__idx_Advection_process_ghosts=0;\
  int Advection::__idx_Advection_doWork=0;\

typedef CBaseT1<ArrayElementT<QuadIndex>, CProxy_Advection> CBase_Advection;

extern void _registerAdvection(void);
#endif
