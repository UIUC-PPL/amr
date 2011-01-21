#ifndef _DECL_Main_H_
#define _DECL_Main_H_
#include "charm++.h"
/* DECLS: readonly CProxy_Main mainProxy;
 */

/* DECLS: readonly int block_height;
 */

/* DECLS: readonly int block_width;
 */

/* DECLS: readonly int array_height;
 */

/* DECLS: readonly int array_width;
 */

/* DECLS: readonly int num_chare_rows;
 */

/* DECLS: readonly int num_chare_cols;
 */

/* DECLS: readonly int nframe;
 */

/* DECLS: readonly double xmin;
 */

/* DECLS: readonly double xmax;
 */

/* DECLS: readonly double ymin;
 */

/* DECLS: readonly double ymax;
 */

/* DECLS: readonly double xctr;
 */

/* DECLS: readonly double yctr;
 */

/* DECLS: readonly double radius;
 */

/* DECLS: readonly int nx;
 */

/* DECLS: readonly int ny;
 */

/* DECLS: readonly double dx;
 */

/* DECLS: readonly double dy;
 */

/* DECLS: readonly double v;
 */

/* DECLS: readonly double ap;
 */

/* DECLS: readonly double an;
 */

/* DECLS: readonly double tmax;
 */

/* DECLS: readonly double dt;
 */

/* DECLS: readonly double cfl;
 */

/* DECLS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
};
 */
 class Main;
 class CkIndex_Main;
 class CProxy_Main;
/* --------------- index object ------------------ */
class CkIndex_Main:public CProxy_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static int __idx_Main_CkArgMsg;
    static int ckNew(CkArgMsg* impl_msg) { return __idx_Main_CkArgMsg; }
    static void _call_Main_CkArgMsg(void* impl_msg,Main* impl_obj);

};
/* --------------- element proxy ------------------ */
class CProxy_Main:public CProxy_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    CProxy_Main(void) {};
    CProxy_Main(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Main(const Chare *c) : CProxy_Chare(c){  }
int ckIsDelegated(void) const {return CProxy_Chare::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_Chare::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_Chare::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_Chare::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_Chare::ckCheck();}
const CkChareID &ckGetChareID(void) const
{ return CProxy_Chare::ckGetChareID(); }
operator const CkChareID &(void) const {return ckGetChareID();}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_Chare::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_Chare::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_Chare::pup(p);
    }
    void ckSetChareID(const CkChareID &c) {
      CProxy_Chare::ckSetChareID(c);
    }
    Main *ckLocal(void) const
     { return (Main *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static CkChareID ckNew(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);
    static void ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int onPE=CK_PE_ANY);
    CProxy_Main(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);

};
PUPmarshall(CProxy_Main)
typedef CBaseT1<Chare, CProxy_Main> CBase_Main;

#include "Advection.decl.h"

extern void _registerMain(void);
extern "C" void CkRegisterMainModule(void);
#endif
