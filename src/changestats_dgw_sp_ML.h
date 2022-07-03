/*  File src/changestats_dgw_sp_ML.h in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_ML_H_
#define _CHANGESTATS_DGW_SP_ML_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm.multi_changestat_multilayer.h"

#define ergm_Check2Path(e11, e12, e21, e22, any_order) ((any_order) ? (((e11)&&(e22)) || ((e12)&&(e21))) : ((e11)&&(e22)))

// FIXME: Optimize
/*! @function
  @abstract Evaluate change in the state of a cross-layer two-path as a result of a particular toggle
  @param tail1,head1,tail2,head2 tails and heads of constituents segments, on layer scale
  @param ll1,ll2 Layer Logics of constituent two-paths
  @param any_order Does which segment is in which layer matter?
  @param c11,c12,c21,c22 cIJ is the precomputed direction of change in the Ith segment's Jth layer value (-1, 0, or +1).
  
  @return -1, 0, or +1
*/


static inline int ergm_c_LayerLogic2Path(Vertex tail1, Vertex head1, Vertex tail2, Vertex head2,
					 StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order,
					 int c11, int c12, int c21, int c22){
  if(!ll1->onwp->directed_flag) any_order = TRUE;
  
  unsigned int e11o = c11? c11==-1 : ML_GETWT(ll1, tail1, head1), e22o = c22? c22==-1 : ML_GETWT(ll2, tail2, head2), e12o=0, e21o=0;
  unsigned int e11n = c11? c11==+1 : e11o, e22n = c22? c22==+1 : e22o, e12n=0, e21n=0;
  
  if(any_order){
    e12o = c12? c12==-1 : ML_GETWT(ll2, tail1, head1);
    e21o = c21? c21==-1 : ML_GETWT(ll1, tail2, head2);
    e12n = c12? c12==+1 : e12o;
    e21n = c21? c21==+1 : e21o;
  }

  return ergm_Check2Path(e11n, e12n, e21n, e22n, any_order)
    - ergm_Check2Path(e11o, e12o, e21o, e22o, any_order);
}


static inline unsigned int ergm_LayerLogic2Path(Vertex tail1, Vertex head1, Vertex tail2, Vertex head2, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order){
  if(!ll1->onwp->directed_flag) any_order = TRUE;

  unsigned int e11 = ML_GETWT(ll1, tail1, head1), e22 = ML_GETWT(ll2, tail2, head2), e12, e21;
  if(any_order){
    e12 = ML_GETWT(ll2, tail1, head1);
    e21 = ML_GETWT(ll1, tail2, head2);
  }
  return ergm_Check2Path(e11, e12, e21, e22, any_order);
}

#define SETUP_update_spcache						\
  Vertex t0 = ML_IO_TAIL(ll0, tail), h0 = ML_IO_HEAD(ll0, head);	\
  int l1fc = ergm_LayerLogic2(t0, h0, tail, head, ll1, TRUE);		\
  int l2fc = ergm_LayerLogic2(t0, h0, tail, head, ll2, TRUE);		\
  int l1rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll1, TRUE) : 0; \
  int l2rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll2, TRUE) : 0; \
  int l3fc = 0, l3rc = 0;


#define CALC_with_dirs(subroutine)					\
  if(l1fc || l2fc || l3fc){						\
    int l1c = l1fc, l2c = l2fc, l3c = l3fc;				\
    (void) l3c;                                                         \
    Vertex t = t0, h = h0;						\
    subroutine;								\
  }									\
  if(l1rc || l2rc || l3rc){						\
    int l1c = l1rc, l2c = l2rc, l3c = l3rc;				\
    (void) l3c;                                                         \
    Vertex t = h0, h = t0;						\
    subroutine;								\
  }

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5

/* /\*DSP calculation functions*\/ */
/* static inline void dspUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* /\* static inline void dspRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); *\/ */
/* static inline void dspOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */


/* /\*ESP calculation functions*\/ */
/* static inline void espUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* /\* static inline void espRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); *\/ */
/* static inline void espOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */

/* /\*Changescore functions*\/ */
/* C_CHANGESTAT_FN(c_desp_ML); */
/* I_CHANGESTAT_FN(i_dgwesp_ML); */
/* C_CHANGESTAT_FN(c_dgwesp_ML); */

/* /\*Changescore functions*\/ */
/* C_CHANGESTAT_FN(c_ddsp_ML); */
/* I_CHANGESTAT_FN(i_dgwdsp_ML); */
/* C_CHANGESTAT_FN(c_dgwdsp_ML); */

/* /\*Changescore functions*\/ */
/* I_CHANGESTAT_FN(i_dnsp_ML); */
/* C_CHANGESTAT_FN(c_dnsp_ML); */
/* I_CHANGESTAT_FN(i_dgwnsp_ML); */
/* C_CHANGESTAT_FN(c_dgwnsp_ML); */

#endif // _CHANGESTATS_DGW_SP_ML_H_
