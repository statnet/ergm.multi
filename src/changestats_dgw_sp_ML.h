/*  File src/changestats_dgw_sp_ML.h in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_ML_H_
#define _CHANGESTATS_DGW_SP_ML_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm.multi_changestat_multilayer.h"
#include "ergm_dyad_hashmap.h"

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
  int l1fc = ergm_LayerLogic2(t0, h0, tail, head, ll1, LL_DIFF);		\
  int l2fc = ergm_LayerLogic2(t0, h0, tail, head, ll2, LL_DIFF);		\
  int l1rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll1, LL_DIFF) : 0; \
  int l2rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll2, LL_DIFF) : 0; \
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


typedef enum {L2UTP, L2OTP, L2ITP, L2RTP, L2OSP, L2ISP} L2Type;

#define call_subroutine_path(count, subroutine_path)    \
  {int L2 = (L2 ## count);                              \
    {subroutine_path}}

#define call_subroutine_focus(count, subroutine_focus)  \
  {int L2 = (L2 ## count);                              \
    {subroutine_focus}}


#define SETUP_calc_dsp							\
  SETUP_update_spcache;

#define SETUP_calc_esp							\
  SETUP_calc_dsp;							\
  l3fc = ergm_LayerLogic2(t0, h0, tail, head, ll3, LL_DIFF);		\
  l3rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll3, LL_DIFF) : 0;


#define INC_IF_TWOPATH(ij, t1, h1, t2, h2) if(ergm_LayerLogic2Path(t1,h1,t2,h2, ll1, ll2, any_order)) L2 ## ij ++;

#define UPDATE_CS_1(ij, t1, h1, t2, h2, subroutine_path)        \
  {                                                             \
    int c2path = ergm_c_LayerLogic2Path(t1,h1,t2,h2,            \
                                        ll1,ll2, any_order,     \
                                        l1c,l2c,0,0);           \
    call_subroutine_path(ij, subroutine_path);                  \
  }

#define UPDATE_CS_2(ij, t1, h1, t2, h2, subroutine_path)        \
  {                                                             \
    int c2path = ergm_c_LayerLogic2Path(t1,h1,t2,h2,            \
                                        ll1,ll2, any_order,     \
                                        0,0,l1c,l2c);           \
    call_subroutine_path(ij, subroutine_path);                  \
  }


/**************************
 dsp Calculation functions
**************************/

/*
  Changescore for ESPs based on two-paths in undirected graphs i.e. configurations for edge i<->j such that i<->k<->j (where <-> here denotes an undirected edge).

  UTP:
  L2th - count t<->k<->h
  L2tk - for each t<->k neq h: k<->h, count u such that k<->u<->h
  L2hk - for each h<->k neq t: k<->t, count u such that k<->u<->t

  This function will only work properly with undirected graphs, and should only be called in that case.
*/
#define dspUTP_ML_change(subroutine_path, subroutine_focus)             \
  SETUP_calc_dsp;                                                       \
  any_order = TRUE;                                                     \
                                                                        \
  CALC_with_dirs({                                                      \
      /* step through edges of head */                                  \
      ML_EXEC_THROUGH_EDGES(ll0, h,e,u, {                               \
          if (u!=t){                                                    \
            unsigned int L2tu = 0;                                      \
            if(spcache) L2tu = GETUDMUI(t,u,spcache);                    \
            else{                                                       \
              /* step through edges of u */                             \
              ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {                       \
                  /* Confirm that 2-path u - v - t satisfies layer conditions */ \
                  INC_IF_TWOPATH(tu,t,v,v,u);                           \
                });                                                     \
            }                                                           \
            UPDATE_CS_1(tu,t,h,h,u,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      ML_EXEC_THROUGH_EDGES(ll0, t,e,u, {                               \
          if (u!=h){                                                    \
            unsigned int L2uh = 0;                                      \
            if(spcache) L2uh = GETUDMUI(u,h,spcache);                    \
            else{                                                       \
              /* step through edges of u */                             \
              ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {                       \
                  /* Confirm that 2-path u - v - h satisfies layer conditions */ \
                  INC_IF_TWOPATH(uh,u,v,v,h);                           \
                });                                                     \
            }                                                           \
            UPDATE_CS_2(uh,u,t,t,h,subroutine_path);                                    \
          }                                                             \
        });                                                             \
    });


/*
  Changescore for dsps based on outgoing two-paths, i.e. configurations for non-edge i->j such that i->k->j.

  This function should only be used in the directed case
*/
#define dspOTP_ML_change(subroutine_path, subroutine_focus)               \
  SETUP_calc_dsp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
      /* step through outedges of head (i.e., k: t->k)*/                \
      ML_EXEC_THROUGH_OUTEDGES(ll0,h, e, k, {                           \
          if(k!=t){ /*Only use contingent cases*/                       \
            unsigned int L2tk = 0;                                      \
            if(spcache) L2tk = GETDDMUI(t,k,spcache);                    \
            else{                                                       \
              /* step through inedges of k, incl. (h,k) itself */       \
              ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {                    \
                  INC_IF_TWOPATH(tk,t,u,u,k);                           \
                });                                                     \
            }                                                           \
            /*Update the changestat for the t->k edge*/                 \
            UPDATE_CS_1(tk,t,h,h,k,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      /* step through inedges of tail (i.e., k: k->t)*/                 \
      ML_EXEC_THROUGH_INEDGES(ll0,t, e, k, {                            \
          if (k!=h){ /*Only use contingent cases*/                      \
            unsigned int L2kh = 0;                                      \
            if(spcache) L2kh = GETDDMUI(k,h,spcache);                    \
            else{                                                       \
              /* step through outedges of k , incl. (k,tail) itself */  \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {                   \
                  INC_IF_TWOPATH(kh,k,u,u,h);                           \
                });                                                     \
            }                                                           \
            /*Update the changestat for the k->t edge*/                 \
            UPDATE_CS_2(kh,k,t,t,h,subroutine_path);                                    \
          }                                                             \
    });                                                                 \
  });


/*
  Changescore for DSPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.
  IE cyclical shared partners

  ITP:
  L2th - count j->k->i
  L2hk - for each j->k neq i: k->i, count u such that k->u->j
  L2kt - for each k->i neq j: j->k, count u such that i->u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspITP_ML_change(subroutine_path, subroutine_focus)               \
  SETUP_calc_dsp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
      /* step through outedges of head (i.e., k: h->k)*/                \
      ML_EXEC_THROUGH_OUTEDGES(ll0, h, e, k, {                          \
          if((k!=t)){ /*Only use contingent cases*/                     \
            unsigned int L2kt = 0;                                      \
            if(spcache) L2kt = GETDDMUI(t,k,spcache); /* spcache is an OTP cache. */ \
            else{                                                       \
              ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {                    \
                  INC_IF_TWOPATH(kt,t,u,u,k);                           \
                });                                                     \
            }                                                           \
            /*Update the changestat for the h->k edge*/                 \
            UPDATE_CS_1(kt,t,h,h,k,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      /* step through inedges of tail (i.e., k: k->t)*/                 \
      ML_EXEC_THROUGH_INEDGES(ll0,t, e, k, {                            \
          if((k!=h)){ /*Only use contingent cases*/                     \
            unsigned int L2hk = 0;                                      \
            if(spcache) L2hk = GETDDMUI(k,h,spcache); /* spcache is an OTP cache. */ \
            else{                                                       \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {                   \
                  INC_IF_TWOPATH(hk,k,u,u,h);                           \
                });                                                     \
            }                                                           \
            /*Update the changestat for the k->t edge*/                 \
            UPDATE_CS_2(hk,k,t,t,h,subroutine_path);                                    \
          }                                                             \
        });                                                             \
    });


/*
  Changescore for DSPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  OSP:
  L2th - count t->k, h->k
  L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
  L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspOSP_ML_change(subroutine_path, subroutine_focus)       \
  SETUP_calc_dsp;                                               \
  any_order = TRUE;                                             \
                                                                \
  CALC_with_dirs({                                              \
      ML_EXEC_THROUGH_INEDGES(ll0,h, e, k, {                    \
          if(k!=t){                                             \
            unsigned int L2tk = 0;                              \
            if(spcache) L2tk = GETUDMUI(t,k,spcache);            \
            else{                                               \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {           \
                  if (u!=t)                                     \
                    /*Increment if there is an OSP  */          \
                    INC_IF_TWOPATH(tk,t,u,k,u);                 \
                });                                             \
            }                                                   \
            /*Update the changestat for the t->k edge*/         \
            /* twice (one for each direction of t<->k) */       \
            UPDATE_CS_1(tk,t,h,k,h,subroutine_path);                            \
          }                                                     \
        });                                                     \
    });


/*
  Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  ISP:
  L2th - count k->t, k->h
  L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
  L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspISP_ML_change(subroutine_path, subroutine_focus)       \
  SETUP_calc_dsp;                                               \
  any_order = TRUE;                                             \
                                                                \
  CALC_with_dirs({                                              \
      ML_EXEC_THROUGH_OUTEDGES(ll0,t, e, k, {                   \
          if(k!=h){                                             \
            unsigned int L2kh = 0;                              \
            if(spcache) L2kh = GETUDMUI(k,h,spcache);            \
            else{                                               \
              ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {            \
                  if(u!=h)                                      \
                    /*Increment if there is an ISP*/            \
                    INC_IF_TWOPATH(kh,u,k,u,h);                 \
                });                                             \
            }                                                   \
            /*Update the changestat for the k->h edge*/         \
            /* twice (one for each direction of h<->k) */       \
            UPDATE_CS_1(kh,t,h,t,k,subroutine_path);                            \
          }                                                     \
        });                                                     \
    });


/* /\* */
/* Changescore for ESPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j). */

/* RTP: */
/* L2th - count t<->k<->h */
/* L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h */
/* L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h */

/* We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function. */
/* *\/ */
/* #define dspRTP_ML_change(subroutine_path, subroutine_focus) \ */
/*   Edge e, f; */
/*   int j, echange, htedge; */
/*   int L2th,L2tk,L2kt,L2hk,L2kh; /\*Two-path counts for various edges*\/ */
/*   Vertex deg; */
/*   Vertex k, u; */

/*   memset(cs, 0, nd*sizeof(double)); */
/*     L2th=0; */
/*     echange = (IS_OUTEDGE(tail,head) == 0) ? 1 : -1; */
/*     htedge=IS_OUTEDGE(head,tail);  /\*Is there an h->t (reciprocating) edge?*\/ */
/*     /\* step through inedges of tail (k->t: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         /\*Do we have a t<->k<->h TP?  If so, add it to our count.*\/ */
/*         L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)); */
/*         if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /\*Only consider stats that could change*\/ */
/*           L2kt=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->t edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kt + echange == deg) - (L2kt == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of tail (t->k: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /\*Only consider stats that could change*\/ */
/*           L2tk=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the t->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2tk + echange == deg) - (L2tk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through inedges of head (k->h: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /\*Only consider stats that could change*\/ */
/*           L2kh=0; */
/*           /\*Now, count # u such that k<->u<->h (to get k->h's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->h edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kh + echange == deg) - (L2kh == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of head (h->k: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /\*Only consider stats that could change*\/ */
/*           L2hk=0; */
/*           /\*Now, count # u such that k<->u<->h (to get h->k's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the h->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2hk + echange == deg) - (L2hk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\*Finally, update the changestat for the t->h edge*\/ */
/*     for(unsigned int j = 0; j < nd; j++){ */
/*       Vertex deg = dvec[j]; */
/*       cs[j] += echange*(L2th == deg); */
/*     } */
/* } */



/**************************
 ESP Calculation functions
**************************/

/*
  Changescore for ESPs based on two-paths in undirected graphs i.e. configurations for edge i<->j such that i<->k<->j (where <-> here denotes an undirected edge).

  UTP:
  L2th - count t<->k<->h
  L2tk - for each t<->k neq h: k<->h, count u such that k<->u<->h
  L2hk - for each h<->k neq t: k<->t, count u such that k<->u<->t

  This function will only work properly with undirected graphs, and should only be called in that case.
*/
#define espUTP_ML_change(subroutine_path, subroutine_focus) \
  SETUP_calc_esp;                                               \
  any_order = TRUE;                                             \
                                                                \
  CALC_with_dirs({                                              \
      unsigned int L2th = 0;                                    \
      if(spcache) L2th = GETUDMUI(t,h,spcache);                  \
      ML_EXEC_THROUGH_EDGES(ll0, h,e,u, {                       \
          if (ML_IS_UNDIRECTED_EDGE(ll0,u,t) != 0){             \
            if(!spcache && l3c) INC_IF_TWOPATH(th,t,u,u,h);     \
            unsigned int Buh = ML_GETWT(ll3,u,h);               \
            unsigned int Btu = ML_GETWT(ll3,t,u);               \
            unsigned int L2tu = 0;                              \
            unsigned int L2uh = 0;                              \
            if(spcache){                                        \
              L2tu = GETUDMUI(t,u,spcache);                      \
              L2uh = GETUDMUI(u,h,spcache);                      \
            }else{                                              \
              ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {               \
                  if(Buh) INC_IF_TWOPATH(uh,u,v,v,h);           \
                  if(Btu) INC_IF_TWOPATH(tu,t,v,v,u);           \
                });                                             \
            }                                                   \
            if(Buh) UPDATE_CS_2(uh,u,t,t,h,subroutine_path);                    \
            if(Btu) UPDATE_CS_1(tu,t,h,h,u,subroutine_path);                    \
          }                                                     \
        });                                                     \
      if(l3c) call_subroutine_focus(th, subroutine_focus);      \
    });


/*
  Changescore for ESPs based on outgoing two-paths, i.e. configurations for edge i->j such that i->k->j.

  OTP:
  L2th - count i->k->j
  L2tk - for each i->k neq j: j->k, count u such that i->u->k
  L2kh - for each k->j neq i: k->i, count u such that k->u->j

  This function should only be used in the directed case, with espUTP being used in the undirected case.
*/
#define espOTP_ML_change(subroutine_path, subroutine_focus) \
  SETUP_calc_esp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
  if(l3c){                                                              \
    unsigned int L2th = 0;                                              \
    if(spcache) L2th = GETDDMUI(t,h,spcache);                            \
    else{                                                               \
      ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, {                             \
	  if(k!=h)                                                      \
	    INC_IF_TWOPATH(th,t,k,k,h);                                 \
	});                                                             \
    }                                                                   \
    call_subroutine_focus(th, subroutine_focus);                        \
  }                                                                     \
  ML_EXEC_THROUGH_OUTEDGES(ll3,t,e,k, {                                 \
  if((k!=h)&&(ML_IS_OUTEDGE(ll0,h,k))){ /*Only use contingent cases*/   \
    unsigned int L2tk = 0;                                              \
    if(spcache) L2tk = GETDDMUI(t,k,spcache);                            \
    else{                                                               \
      ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {                              \
          if(u!=t)                                                      \
            INC_IF_TWOPATH(tk,t,u,u,k);                                 \
        });                                                             \
    }                                                                   \
    UPDATE_CS_1(tk,t,h,h,k,subroutine_path);                                            \
  }                                                                     \
    });                                                                 \
  /* step through inedges of h (i.e., k: k->h)*/                        \
  ML_EXEC_THROUGH_INEDGES(ll3,h,e,k, {                                  \
      if((k!=t)&&(ML_IS_OUTEDGE(ll0,k,t))){ /*Only use contingent cases*/ \
        unsigned int L2kh = 0;                                          \
	if(spcache) L2kh = GETDDMUI(k,h,spcache);                        \
	else{                                                           \
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {                         \
	      if(u!=h)                                                  \
		INC_IF_TWOPATH(kh,k,u,u,h);                             \
	    });                                                         \
	}                                                               \
	UPDATE_CS_2(kh,k,t,t,h,subroutine_path);                                        \
      }                                                                 \
    });                                                                 \
  });


/*
  Changescore for ESPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.

  ITP:
  L2th - count j->k->i
  L2hk - for each j->k neq i: k->i, count u such that k->u->j
  L2kt - for each k->i neq j: j->k, count u such that i->u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espITP_ML_change(subroutine_path, subroutine_focus)               \
  SETUP_calc_esp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
      if(l3c){                                                          \
        unsigned int L2th = 0;                                          \
        if(spcache) L2th = GETDDMUI(h,t,spcache); /* spcache is an OTP cache. */ \
        else{                                                           \
          ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, {                          \
              if(k!=h)                                                  \
                INC_IF_TWOPATH(th,h,k,k,t);                             \
            });                                                         \
        }                                                               \
        call_subroutine_focus(th, subroutine_focus);                    \
      }                                                                 \
                                                                        \
      ML_EXEC_THROUGH_OUTEDGES(ll3,h,e,k, {                             \
          if((k!=t)&&(ML_IS_OUTEDGE(ll0,k,t))){ /*Only use contingent cases*/ \
            /*We have a h->k->t two-path, so add it to our count.*/     \
            unsigned int L2hk = 0;                                      \
            if(spcache) L2hk = GETDDMUI(k,h,spcache); /* spcache is an OTP cache. */ \
            else{                                                       \
              /*Now, count # u such that k->u->h (so that we know k's ESP value)*/ \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {                     \
                  if(u!=h)                                              \
                    INC_IF_TWOPATH(hk,k,u,u,h);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the h->k edge*/                 \
            UPDATE_CS_2(hk,k,t,t,h,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      /* step through inedges of t (i.e., k: k->t)*/                    \
      ML_EXEC_THROUGH_INEDGES(ll3,t,e,k, {                              \
          if((k!=h)&&(ML_IS_OUTEDGE(ll0,h,k))){ /*Only use contingent cases*/ \
            unsigned int L2kt = 0;                                      \
            if(spcache) L2kt = GETDDMUI(t,k,spcache); /* spcache is an OTP cache. */ \
            else{                                                       \
              /*Now, count # u such that t->u->k (so that we know k's ESP value)*/ \
              ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {                      \
                  if(u!=t)                                              \
                    INC_IF_TWOPATH(kt,t,u,u,k);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the k->t edge*/                 \
            UPDATE_CS_1(kt,t,h,h,k,subroutine_path);                                    \
          }                                                             \
        });                                                             \
    });


/*
  Changescore for ESPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  OSP:
  L2th - count t->k, h->k
  L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
  L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espOSP_ML_change(subroutine_path, subroutine_focus)               \
  SETUP_calc_esp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
      if(l3c){                                                          \
        unsigned int L2th = 0;                                          \
        if(spcache) L2th = GETUDMUI(t,h,spcache);                        \
        else{                                                           \
          ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, {                         \
              if(k!=h)                                                  \
                INC_IF_TWOPATH(th,t,k,h,k);                             \
            });                                                         \
        }                                                               \
        call_subroutine_focus(th, subroutine_focus);                    \
      }                                                                 \
                                                                        \
      /* step through outedges of t (i.e., k: t->k, k->h, k!=h)*/       \
      ML_EXEC_THROUGH_OUTEDGES(ll3,t,e,k, {                             \
          if(k!=h && ML_IS_OUTEDGE(ll0,k,h)){ /*Only consider stats that could change*/ \
            unsigned int L2tk = 0;                                      \
            if(spcache) L2tk = GETUDMUI(t,k,spcache);                    \
            else{                                                       \
              /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/ \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {                     \
                  if(u!=t)                                              \
                    /*Increment if there is an OSP*/                    \
                    INC_IF_TWOPATH(tk,t,u,k,u);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the t->k edge*/                 \
            UPDATE_CS_1(tk,t,h,k,h,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      /* step through inedges of t (i.e., k: k->t, k->h, k!=h)*/        \
      ML_EXEC_THROUGH_INEDGES(ll3,t,e,k, {                              \
          if(k!=h && ML_IS_OUTEDGE(ll0,k,h)){ /*Only stats that could change*/ \
            unsigned int L2kt=0;                                        \
            if(spcache) L2kt = GETUDMUI(k,t,spcache);                    \
            else{                                                       \
              /*Now, count # u such that t->u,k->u (to get k->t's ESP value)*/ \
              ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {                     \
                  if(u!=t)                                              \
                    /*Increment if there is an OSP*/                    \
                    INC_IF_TWOPATH(kt,t,u,k,u);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the k->t edge*/                 \
            UPDATE_CS_1(kt,t,h,k,h,subroutine_path);                                    \
          }                                                             \
        });                                                             \
    });


/*
  Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  ISP:
  L2th - count k->t, k->h
  L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
  L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espISP_ML_change(subroutine_path, subroutine_focus)               \
  SETUP_calc_esp;                                                       \
                                                                        \
  CALC_with_dirs({                                                      \
      if(l3c){                                                          \
        unsigned int L2th = 0;                                          \
        if(spcache) L2th = GETUDMUI(t,h,spcache);                        \
        else{                                                           \
          ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, {                          \
              if(k!=h)                                                  \
                INC_IF_TWOPATH(th,k,t,k,h);                             \
            });                                                         \
        }                                                               \
        call_subroutine_focus(th, subroutine_focus);                    \
      }                                                                 \
                                                                        \
      /* step through inedges of h (i.e., k: k->h, t->k, k!=t)*/        \
      ML_EXEC_THROUGH_INEDGES(ll3,h,e,k, {                              \
          if(k!=t && ML_IS_OUTEDGE(ll0,t,k)){                           \
            unsigned int L2kh = 0;                                      \
            if(spcache) L2kh = GETUDMUI(k,h,spcache);                    \
            else{                                                       \
              /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/ \
              ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {                      \
                  if(u!=h)                                              \
                    /*Increment if there is an ISP*/                    \
                    INC_IF_TWOPATH(kh,u,k,u,h);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the k->h edge*/                 \
            UPDATE_CS_1(kh,t,h,t,k,subroutine_path);                                    \
          }                                                             \
        });                                                             \
      /* step through outedges of h (i.e., k: h->k, t->k, k!=t)*/       \
      ML_EXEC_THROUGH_OUTEDGES(ll3,h,e,k, {                             \
          if(k!=t && ML_IS_OUTEDGE(ll0,t,k)){ /*Only stats that could change*/ \
            unsigned int L2hk = 0;                                      \
            if(spcache) L2hk = GETUDMUI(h,k,spcache);                    \
            else{                                                       \
              /*Now, count # u such that u->h,u->k (to get k->h's ESP value)*/ \
              ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {                      \
                  if(u!=h)                                              \
                    /*Increment if there is an ISP*/                    \
                    INC_IF_TWOPATH(hk,u,k,u,h);                         \
                });                                                     \
            }                                                           \
            /*Update the changestat for the h->k edge*/                 \
            UPDATE_CS_1(hk,t,h,t,k,subroutine_path);                                    \
          }                                                             \
        });                                                             \
    });


/* /\* */
/* Changescore for ESPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j). */

/* RTP: */
/* L2th - count t<->k<->h */
/* L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h */
/* L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h */

/* We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function. */
/* *\/ */
/* #define espRTP_ML_change(subroutine_path, subroutine_focus) \ */
/*   Edge e, f; */
/*   int j, echange, htedge; */
/*   int L2th,L2tk,L2kt,L2hk,L2kh; /\*Two-path counts for various edges*\/ */
/*   Vertex deg; */
/*   Vertex k, u; */

/*   memset(cs, 0, nd*sizeof(double)); */

/*     L2th=0; */
/*     echange = (IS_OUTEDGE(t,h) == 0) ? 1 : -1; */
/*     htedge=IS_OUTEDGE(h,t);  /\*Is there an h->t (reciprocating) edge?*\/ */
/*     /\* step through inedges of t (k->t: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         /\*Do we have a t<->k<->h TP?  If so, add it to our count.*\/ */
/*         L2th+=(IS_OUTEDGE(t,k)&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)); */
/*         if(htedge&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)){ /\*Only consider stats that could change*\/ */
/*           L2kt=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2kt+=(IS_OUTEDGE(u,t)&&IS_OUTEDGE(t,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->t edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kt + echange == deg) - (L2kt == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of t (t->k: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         if(htedge&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)){ /\*Only consider stats that could change*\/ */
/*           L2tk=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2tk+=(IS_OUTEDGE(u,t)&&IS_OUTEDGE(t,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the t->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2tk + echange == deg) - (L2tk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through inedges of h (k->h: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(t,k)&&IS_OUTEDGE(k,t)){ /\*Only consider stats that could change*\/ */
/*           L2kh=0; */
/*           /\*Now, count # u such that k<->u<->h (to get k->h's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2kh+=(IS_OUTEDGE(u,h)&&IS_OUTEDGE(h,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->h edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kh + echange == deg) - (L2kh == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of h (h->k: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(t,k)&&IS_OUTEDGE(k,t)){ /\*Only consider stats that could change*\/ */
/*           L2hk=0; */
/*           /\*Now, count # u such that k<->u<->h (to get h->k's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2hk+=(IS_OUTEDGE(u,h)&&IS_OUTEDGE(h,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the h->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2hk + echange == deg) - (L2hk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\*Finally, update the changestat for the t->h edge*\/ */
/*     for(unsigned int j = 0; j < nd; j++){ */
/*       Vertex deg = dvec[j]; */
/*       cs[j] += echange*(L2th == deg); */
/*     } */
/* } */


#endif // _CHANGESTATS_DGW_SP_ML_H_
