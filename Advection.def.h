/* DEFS: array Advection: ArrayElement{
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
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Advection::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Advection(CkMigrateMessage* impl_msg);
 */

/* DEFS: Advection(void);
 */
void CProxyElement_Advection::insert(int onPE)
{ 
  void *impl_msg = CkAllocSysMsg();
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Advection::__idx_Advection_void,onPE);
}

/* DEFS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */
void CProxyElement_Advection::insert(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2, int onPE, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_2;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_2;
  }
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Advection::__idx_Advection_marshall2,onPE);
}

/* DEFS: void refine(void);
 */
void CProxyElement_Advection::refine(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_refine_void,0);
}

/* DEFS: void inform_nbr_of_refinement(void);
 */
void CProxyElement_Advection::inform_nbr_of_refinement(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_inform_nbr_of_refinement_void,0);
}

/* DEFS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
void CProxyElement_Advection::receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int iter, int dir, int width, const double *u
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_u, impl_cnt_u;
  impl_off_u=impl_off=CK_ALIGN(impl_off,sizeof(double));
  impl_off+=(impl_cnt_u=sizeof(double)*(width));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_u,u,impl_cnt_u);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_receiveGhosts_marshall5,0);
}

/* DEFS: void begin_iteration(void);
 */
void CProxyElement_Advection::begin_iteration(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_begin_iteration_void,0);
}

/* DEFS: void doStep(void);
 */
void CProxyElement_Advection::doStep(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_doStep_void,0);
}

/* DEFS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
void CProxyElement_Advection::requestNextFrame(liveVizRequestMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_requestNextFrame_liveVizRequestMsg,0);
}

/* DEFS: void setReal(void);
 */
void CProxyElement_Advection::setReal(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_setReal_void,0);
}

/* DEFS: void setDummy(void);
 */
void CProxyElement_Advection::setDummy(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_setDummy_void,0);
}

/* DEFS: Advection(CkMigrateMessage* impl_msg);
 */
 int CkIndex_Advection::__idx_Advection_CkMigrateMessage=0;
void CkIndex_Advection::_call_Advection_CkMigrateMessage(void* impl_msg,Advection * impl_obj)
{
  new (impl_obj) Advection((CkMigrateMessage*)impl_msg);
}

/* DEFS: Advection(void);
 */
CkArrayID CProxy_Advection::ckNew(const CkArrayOptions &opts)
{ 
  void *impl_msg = CkAllocSysMsg();
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_Advection::__idx_Advection_void,opts);
}
 int CkIndex_Advection::__idx_Advection_void=0;
void CkIndex_Advection::_call_Advection_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) Advection();
}

/* DEFS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */
CkArrayID CProxy_Advection::ckNew(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_2;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(bool &)impl_noname_2;
  }
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_Advection::__idx_Advection_marshall2,opts);
}
 int CkIndex_Advection::__idx_Advection_marshall2=0;
void CkIndex_Advection::_call_Advection_marshall2(void* impl_msg,Advection * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  bool impl_noname_0; implP|impl_noname_0;
  bool impl_noname_1; implP|impl_noname_1;
  bool impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) Advection(impl_noname_0, impl_noname_1, impl_noname_2);
}
int CkIndex_Advection::_callmarshall_Advection_marshall2(char* impl_buf,Advection * impl_obj) {
  /*Unmarshall pup'd fields: const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  bool impl_noname_0; implP|impl_noname_0;
  bool impl_noname_1; implP|impl_noname_1;
  bool impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) Advection(impl_noname_0, impl_noname_1, impl_noname_2);
  return implP.size();
}
void CkIndex_Advection::_marshallmessagepup_Advection_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  bool impl_noname_0; implP|impl_noname_0;
  bool impl_noname_1; implP|impl_noname_1;
  bool impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
}

/* DEFS: void refine(void);
 */
void CProxy_Advection::refine(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_refine_void,0);
}
 int CkIndex_Advection::__idx_refine_void=0;
void CkIndex_Advection::_call_refine_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->refine();
}

/* DEFS: void inform_nbr_of_refinement(void);
 */
void CProxy_Advection::inform_nbr_of_refinement(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_inform_nbr_of_refinement_void,0);
}
 int CkIndex_Advection::__idx_inform_nbr_of_refinement_void=0;
void CkIndex_Advection::_call_inform_nbr_of_refinement_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->inform_nbr_of_refinement();
}

/* DEFS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
void CProxy_Advection::receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int iter, int dir, int width, const double *u
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_u, impl_cnt_u;
  impl_off_u=impl_off=CK_ALIGN(impl_off,sizeof(double));
  impl_off+=(impl_cnt_u=sizeof(double)*(width));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_u,u,impl_cnt_u);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_receiveGhosts_marshall5,0);
}
 int CkIndex_Advection::__idx_receiveGhosts_marshall5=0;
void CkIndex_Advection::_call_receiveGhosts_marshall5(void* impl_msg,Advection * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int iter, int dir, int width, const double *u*/
  PUP::fromMem implP(impl_buf);
  int iter; implP|iter;
  int dir; implP|dir;
  int width; implP|width;
  int impl_off_u, impl_cnt_u; 
  implP|impl_off_u;
  implP|impl_cnt_u;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  double *u=(double *)(impl_buf+impl_off_u);
  impl_obj->receiveGhosts(iter, dir, width, u);
}
int CkIndex_Advection::_callmarshall_receiveGhosts_marshall5(char* impl_buf,Advection * impl_obj) {
  /*Unmarshall pup'd fields: int iter, int dir, int width, const double *u*/
  PUP::fromMem implP(impl_buf);
  int iter; implP|iter;
  int dir; implP|dir;
  int width; implP|width;
  int impl_off_u, impl_cnt_u; 
  implP|impl_off_u;
  implP|impl_cnt_u;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  double *u=(double *)(impl_buf+impl_off_u);
  impl_obj->receiveGhosts(iter, dir, width, u);
  return implP.size();
}
void CkIndex_Advection::_marshallmessagepup_receiveGhosts_marshall5(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int iter, int dir, int width, const double *u*/
  PUP::fromMem implP(impl_buf);
  int iter; implP|iter;
  int dir; implP|dir;
  int width; implP|width;
  int impl_off_u, impl_cnt_u; 
  implP|impl_off_u;
  implP|impl_cnt_u;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  double *u=(double *)(impl_buf+impl_off_u);
  if (implDestP.hasComments()) implDestP.comment("iter");
  implDestP|iter;
  if (implDestP.hasComments()) implDestP.comment("dir");
  implDestP|dir;
  if (implDestP.hasComments()) implDestP.comment("width");
  implDestP|width;
  if (implDestP.hasComments()) implDestP.comment("u");
  implDestP.synchronize(PUP::sync_begin_array);
  { for (int impl_i=0;impl_i*(sizeof(*u))<impl_cnt_u;impl_i++) { 
      implDestP.synchronize(PUP::sync_item);
      implDestP|u[impl_i];
  } } 
  implDestP.synchronize(PUP::sync_end_array);
}

/* DEFS: void begin_iteration(void);
 */
void CProxy_Advection::begin_iteration(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_begin_iteration_void,0);
}
 int CkIndex_Advection::__idx_begin_iteration_void=0;
void CkIndex_Advection::_call_begin_iteration_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->begin_iteration();
}

/* DEFS: void doStep(void);
 */
void CProxy_Advection::doStep(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_doStep_void,0);
}
 int CkIndex_Advection::__idx_doStep_void=0;
void CkIndex_Advection::_call_doStep_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->doStep();
}

/* DEFS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
void CProxy_Advection::requestNextFrame(liveVizRequestMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_requestNextFrame_liveVizRequestMsg,0);
}
 int CkIndex_Advection::__idx_requestNextFrame_liveVizRequestMsg=0;
void CkIndex_Advection::_call_requestNextFrame_liveVizRequestMsg(void* impl_msg,Advection * impl_obj)
{
  impl_obj->requestNextFrame((liveVizRequestMsg*)impl_msg);
}

/* DEFS: void setReal(void);
 */
void CProxy_Advection::setReal(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_setReal_void,0);
}
 int CkIndex_Advection::__idx_setReal_void=0;
void CkIndex_Advection::_call_setReal_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->setReal();
}

/* DEFS: void setDummy(void);
 */
void CProxy_Advection::setDummy(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Advection::__idx_setDummy_void,0);
}
 int CkIndex_Advection::__idx_setDummy_void=0;
void CkIndex_Advection::_call_setDummy_void(void* impl_msg,Advection * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->setDummy();
}

/* DEFS: Advection(CkMigrateMessage* impl_msg);
 */

/* DEFS: Advection(void);
 */

/* DEFS: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
 */

/* DEFS: void refine(void);
 */
void CProxySection_Advection::refine(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_refine_void,0);
}

/* DEFS: void inform_nbr_of_refinement(void);
 */
void CProxySection_Advection::inform_nbr_of_refinement(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_inform_nbr_of_refinement_void,0);
}

/* DEFS: void receiveGhosts(int iter, int dir, int width, const double *u);
 */
void CProxySection_Advection::receiveGhosts(int iter, int dir, int width, const double *u, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int iter, int dir, int width, const double *u
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_u, impl_cnt_u;
  impl_off_u=impl_off=CK_ALIGN(impl_off,sizeof(double));
  impl_off+=(impl_cnt_u=sizeof(double)*(width));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|iter;
    implP|dir;
    implP|width;
    implP|impl_off_u;
    implP|impl_cnt_u;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_u,u,impl_cnt_u);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_receiveGhosts_marshall5,0);
}

/* DEFS: void begin_iteration(void);
 */
void CProxySection_Advection::begin_iteration(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_begin_iteration_void,0);
}

/* DEFS: void doStep(void);
 */
void CProxySection_Advection::doStep(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_doStep_void,0);
}

/* DEFS: void requestNextFrame(liveVizRequestMsg* impl_msg);
 */
void CProxySection_Advection::requestNextFrame(liveVizRequestMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_requestNextFrame_liveVizRequestMsg,0);
}

/* DEFS: void setReal(void);
 */
void CProxySection_Advection::setReal(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_setReal_void,0);
}

/* DEFS: void setDummy(void);
 */
void CProxySection_Advection::setDummy(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Advection::__idx_setDummy_void,0);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Advection::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
// REG: Advection(CkMigrateMessage* impl_msg);
  __idx_Advection_CkMigrateMessage = CkRegisterEp("Advection(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_Advection_CkMigrateMessage, 0, __idx, 0);
  CkRegisterMigCtor(__idx, __idx_Advection_CkMigrateMessage);

// REG: Advection(void);
  __idx_Advection_void = CkRegisterEp("Advection(void)",
     (CkCallFnPtr)_call_Advection_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_Advection_void);

// REG: Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2);
  __idx_Advection_marshall2 = CkRegisterEp("Advection(const bool &impl_noname_0, const bool &impl_noname_1, const bool &impl_noname_2)",
     (CkCallFnPtr)_call_Advection_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_Advection_marshall2,(CkMarshallUnpackFn)_callmarshall_Advection_marshall2);
  CkRegisterMessagePupFn(__idx_Advection_marshall2,(CkMessagePupFn)_marshallmessagepup_Advection_marshall2);

// REG: void refine(void);
  __idx_refine_void = CkRegisterEp("refine(void)",
     (CkCallFnPtr)_call_refine_void, 0, __idx, 0);

// REG: void inform_nbr_of_refinement(void);
  __idx_inform_nbr_of_refinement_void = CkRegisterEp("inform_nbr_of_refinement(void)",
     (CkCallFnPtr)_call_inform_nbr_of_refinement_void, 0, __idx, 0);

// REG: void receiveGhosts(int iter, int dir, int width, const double *u);
  __idx_receiveGhosts_marshall5 = CkRegisterEp("receiveGhosts(int iter, int dir, int width, const double *u)",
     (CkCallFnPtr)_call_receiveGhosts_marshall5, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_receiveGhosts_marshall5,(CkMarshallUnpackFn)_callmarshall_receiveGhosts_marshall5);
  CkRegisterMessagePupFn(__idx_receiveGhosts_marshall5,(CkMessagePupFn)_marshallmessagepup_receiveGhosts_marshall5);

// REG: void begin_iteration(void);
  __idx_begin_iteration_void = CkRegisterEp("begin_iteration(void)",
     (CkCallFnPtr)_call_begin_iteration_void, 0, __idx, 0);

// REG: void doStep(void);
  __idx_doStep_void = CkRegisterEp("doStep(void)",
     (CkCallFnPtr)_call_doStep_void, 0, __idx, 0);

// REG: void requestNextFrame(liveVizRequestMsg* impl_msg);
  __idx_requestNextFrame_liveVizRequestMsg = CkRegisterEp("requestNextFrame(liveVizRequestMsg* impl_msg)",
     (CkCallFnPtr)_call_requestNextFrame_liveVizRequestMsg, CMessage_liveVizRequestMsg::__idx, __idx, 0);

// REG: void setReal(void);
  __idx_setReal_void = CkRegisterEp("setReal(void)",
     (CkCallFnPtr)_call_setReal_void, 0, __idx, 0);

// REG: void setDummy(void);
  __idx_setDummy_void = CkRegisterEp("setDummy(void)",
     (CkCallFnPtr)_call_setDummy_void, 0, __idx, 0);

  Advection::__sdag_register(); 
}
#endif

Advection_SDAG_CODE_DEF


#ifndef CK_TEMPLATES_ONLY
void _registerAdvection(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: array Advection: ArrayElement{
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
  CkIndex_Advection::__register("Advection", sizeof(Advection));

}
#endif
