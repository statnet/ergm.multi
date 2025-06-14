/*  File inst/include/ergm.multi_changestat_multilayer.h in package ergm.multi,
 *  part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_MULTI_CHANGESTAT_MULTILAYER_H_
#define _ERGM_MULTI_CHANGESTAT_MULTILAYER_H_

#include "ergm_changestat.h"
#include "inc/ergm.multi_changestat_multilayer_common.do_not_include_directly.h"

/* layer-aware macros eponymous to ergm_changestat.h

   The first argument (ll) is a StoreLayerLogic struct defined below.
 */
#define ML_IS_OUTEDGE(ll, a,b) (EdgetreeSearch((a),(b),(ll)->onwp->outedges)!=0?1:0)
#define ML_IS_INEDGE(ll, a,b) (EdgetreeSearch((a),(b),(ll)->onwp->inedges)!=0?1:0)
#define ML_IS_UNDIRECTED_EDGE(ll, a,b) ML_IS_OUTEDGE((ll), MIN(a,b), MAX(a,b))
#define ML_MIN_OUTEDGE(ll, a) (EdgetreeMinimum((ll)->onwp->outedges, (a)))
#define ML_MIN_INEDGE(ll, a) (EdgetreeMinimum((ll)->onwp->inedges, (a)))
#define ML_NEXT_OUTEDGE(ll, e) (EdgetreeSuccessor((ll)->onwp->outedges,(e)))
#define ML_NEXT_INEDGE(ll, e) (EdgetreeSuccessor((ll)->onwp->inedges,(e)))
#define ML_NEXT_OUTEDGE_PRE(ll, e) (EdgetreePreSuccessor((ll)->onwp->outedges,(e)))
#define ML_NEXT_INEDGE_PRE(ll, e) (EdgetreePreSuccessor((ll)->onwp->inedges,(e)))
#define ML_STEP_THROUGH_OUTEDGES(ll, a,e,v) for((e)=ML_MIN_OUTEDGE((ll), a);((v)=ML_OUTVAL((ll), e))!=0;(e)=ML_NEXT_OUTEDGE((ll), e))
#define ML_STEP_THROUGH_INEDGES(ll, a,e,v) for((e)=ML_MIN_INEDGE((ll), a);((v)=ML_INVAL((ll), e))!=0;(e)=ML_NEXT_INEDGE((ll), e))
#define ML_STEP_THROUGH_OUTEDGES_PRE(ll, a,e,v) for((e)=(a);((v)=ML_OUTVAL((ll), e))!=0;(e)=ML_NEXT_OUTEDGE_PRE((ll), e))
#define ML_STEP_THROUGH_INEDGES_PRE(ll, a,e,v) for((e)=(a);((v)=ML_INVAL((ll), e))!=0;(e)=ML_NEXT_INEDGE_PRE((ll), e))
#define ML_STEP_THROUGH_OUTEDGES_DECL(ll, a,e,v) Vertex v; for(Edge e=ML_MIN_OUTEDGE((ll), a);((v)=ML_OUTVAL((ll), e))!=0;e=ML_NEXT_OUTEDGE((ll), e))
#define ML_STEP_THROUGH_INEDGES_DECL(ll, a,e,v) Vertex v; for(Edge e=ML_MIN_INEDGE((ll), a);((v)=ML_INVAL((ll), e))!=0;e=ML_NEXT_INEDGE((ll), e))
#define ML_STEP_THROUGH_OUTEDGES_PRE_DECL(ll, a,e,v) Vertex v; for(Edge e=(a);((v)=ML_OUTVAL((ll), e))!=0;e=ML_NEXT_OUTEDGE_PRE((ll), e))
#define ML_STEP_THROUGH_INEDGES_PRE_DECL(ll, a,e,v) Vertex v; for(Edge e=(a);((v)=ML_INVAL((ll), e))!=0;e=ML_NEXT_INEDGE_PRE((ll), e))
#define ML_EXEC_THROUGH_OUTEDGES(ll, a,e,v,subroutine) {if(ML_DIRECTED((ll))){ ML_EXEC_THROUGH_FOUTEDGES((ll), a,e,v, subroutine) } else { ML_EXEC_THROUGH_EDGES((ll), a,e,v,subroutine) }}
#define ML_EXEC_THROUGH_INEDGES(ll, a,e,v,subroutine) {if(ML_DIRECTED((ll))){ ML_EXEC_THROUGH_FINEDGES((ll), a,e,v, subroutine) } else { ML_EXEC_THROUGH_EDGES((ll), a,e,v,subroutine) }}
#define ML_EXEC_THROUGH_EDGES(ll, a,e,v,subroutine) { ML_EXEC_THROUGH_FOUTEDGES((ll), a,e,v, subroutine) ML_EXEC_THROUGH_FINEDGES((ll), a,e,v, subroutine) }
#define ML_EXEC_THROUGH_OUTEDGES_PRE(ll, a,e,v,subroutine) {if(ML_DIRECTED((ll))){ ML_EXEC_THROUGH_FOUTEDGES_PRE((ll), a,e,v, subroutine) } else { ML_EXEC_THROUGH_EDGES_PRE((ll), a,e,v,subroutine) }}
#define ML_EXEC_THROUGH_INEDGES_PRE(ll, a,e,v,subroutine) {if(ML_DIRECTED((ll))){ ML_EXEC_THROUGH_FINEDGES_PRE((ll), a,e,v, subroutine) } else { ML_EXEC_THROUGH_EDGES_PRE((ll), a,e,v,subroutine) }}
#define ML_EXEC_THROUGH_EDGES_PRE(ll, a,e,v,subroutine) { ML_EXEC_THROUGH_FOUTEDGES_PRE((ll), a,e,v, subroutine) ML_EXEC_THROUGH_FINEDGES_PRE((ll), a,e,v, subroutine) }
#define ML_EXEC_THROUGH_FOUTEDGES(ll, a,e,v,subroutine) {ML_STEP_THROUGH_OUTEDGES_DECL((ll), a,e,v) {subroutine}}
#define ML_EXEC_THROUGH_FINEDGES(ll, a,e,v,subroutine) {ML_STEP_THROUGH_INEDGES_DECL((ll), a,e,v) {subroutine}}
#define ML_EXEC_THROUGH_NET_EDGES(ll, a,b,e,subroutine) {for(Vertex a=1; a <= (ll)->onwp->nnodes; a++){ML_EXEC_THROUGH_FOUTEDGES((ll), a, e, b, {subroutine});}}
#define ML_EXEC_THROUGH_FOUTEDGES_PRE(ll, a,e,v,subroutine) ML_STEP_THROUGH_OUTEDGES_PRE_DECL((ll), a,e,v) {subroutine}
#define ML_EXEC_THROUGH_FINEDGES_PRE(ll, a,e,v,subroutine) ML_STEP_THROUGH_INEDGES_PRE_DECL((ll), a,e,v) {subroutine}
#define ML_EXEC_THROUGH_NET_EDGES_PRE(ll, a,b,e,subroutine) {for(Vertex a=1; a <= (ll)->onwp->nnodes; a++){ML_EXEC_THROUGH_FOUTEDGES_PRE((ll), a, e, b, {subroutine});}}
#define ML_TOGGLE(ll, a,b) (ToggleEdge((a),(b),(ll)->onwp))
#define ML_TOGGLE_DISCORD(ll, a,b) (ToggleEdge((a),(b),(ll)->onwp+1))
#define ML_GETWT(ll, a,b) (GetEdge(a,b,(ll)->onwp))
#define ML_SETWT(ll, a,b,w) (SetEdge(a,b,w,(ll)->onwp))

#define ML_STOP (-INT_MAX)

typedef struct {
  unsigned int nl;
  Network *inwp, *onwp;
  Vertex *lid;
  Vertex *lmap;
  int *symm;
  Rboolean need_ht;
  int *commands;
  int *stack0;
  int *stack1;
} StoreLayerLogic;

/* The vertex index on the input (i.e., combined) network
   corresponding to the vertex on the output (i.e., logical layer)
   network. */
#define ML_OI_TAIL(ll, l, t) ((Vertex) ((ll)->inwp->bipartite? (t) + ((l)-1)*(ll)->onwp->bipartite : (t) + ((l)-1)*(ll)->onwp->nnodes))
#define ML_OI_HEAD(ll, l, h) ((Vertex) ((ll)->inwp->bipartite? (h) + (ll)->inwp->bipartite - (ll)->onwp->bipartite + ((l)-1)*((ll)->onwp->nnodes-(ll)->onwp->bipartite) : (h) + ((l)-1)*(ll)->onwp->nnodes))

/* The vertex index and layer index on the output (i.e., logical
   layer) network corresponding to the vertex on the input (i.e.,
   combined) network. */
#define ML_IO_TAIL(ll, t) ((Vertex) ((ll)->lmap[t]))
#define ML_IO_HEAD(ll, h) ((Vertex) (ll)->lmap[h])
#define ML_LID_TAIL(ll, t) ((Vertex) ((ll)->lid[t]))
#define ML_LID_HEAD(ll, h) ((Vertex) ((ll)->lid[h]))

/* Shorthand for getting and setting the specified edge on the input
   network corresponding to the specified edge on the output
   network. */
#define ML_IGETWT(ll, l,a,b) (GetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b), ll->inwp))
#define ML_ISETWT(ll, l,a,b,w) (SetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b),w,(ll)->inwp))

/*
  Layer evaluation

  The three functions below are used to evaluate logical layers.

  ergm_LayerLogic2() is the workhorse function that takes the focus
    dyad (on the logical layer) and the toggled dyad (on the combined
    network) and returns an int whose value depends on the task:

    LL_ASIS: Just evaluate as is; ttail and thead are ignored.

    LL_DIFF: Difference between a hypothetical toggle of (ttail,thead)
      and its current state.
   
    LL_ENCODE: A binary encoding of the pre- and post-toggle state:
      asis*1 + toggled*2, which can be extracted using a mask
      (output&1 and (output&2)!=0).

    LL_POST: Post-toggle value.

  ergm_LayerLogic() is a shortcut when the same dyad is being looked
    up as being toggled.

  ergm_LayerLogic_affects() has two functions:

    1) If tasked with LL_ENCODE, it encodes the pre- and post-toggle
      state for both the dyad on the logical layer corresponding to
      the toggled dyad on the combined network *and* its opposite
      (head,tail), since some layer logic operations can reference the
      reverse dyad.

    2) Otherwise, it saves the list of dyads on the logical layer
      affected by the toggle and returns their count.
 */

typedef enum{LL_ASIS, LL_DIFF, LL_ENCODE, LL_POST} LayerLogicTask;

/*
  Layer logic language:

  com > 0: reference to layer; look up dyad value and push

  com == 0: numeric literal; push value of next command

  see lookup table in R/InitErgmTerm.multilayer.R
    test_eval.LayerLogic() for com < 0
*/

#define ergm_UNOP(op)				\
  {						\
    int x0 = *(stack0--);			\
    *(++stack0) = (op x0);			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      *(++stack1) = (op x1);			\
    }						\
    break;}

#define ergm_UNFUN(fun)				\
  {						\
    int x0 = *(stack0--);			\
    *(++stack0) = fun(x0);			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      *(++stack1) = fun(x1);			\
    }						\
    break;}


#define ergm_BINOP(op)				\
  {						\
    int x0 = *(stack0--);			\
    int y0 = *(stack0--);			\
    *(++stack0) = (x0 op y0);			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      int y1 = *(stack1--);			\
      *(++stack1) = (x1 op y1);			\
    }						\
    break;}

#define ergm_BINFUN(fun)			\
  {						\
    int x0 = *(stack0--);			\
    int y0 = *(stack0--);			\
    *(++stack0) = fun(x0, y0);			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      int y1 = *(stack1--);			\
      *(++stack1) = fun(x1, y1);		\
    }						\
    break;}

#define ergm_LUNOP(op)				\
  {						\
    int x0 = *(stack0--);			\
    *(++stack0) = (op (x0!=0));			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      *(++stack1) = (op (x1!=0));		\
    }						\
    break;}

#define ergm_LUNFUN(fun)			\
  {						\
    int x0 = *(stack0--);			\
    *(++stack0) = fun(x0!=0);			\
    if(stack1){					\
      int x1 = *(stack1--);			\
      *(++stack1) = fun(x1!=0);			\
    }						\
    break;}


#define ergm_LBINOP(op)				\
  {						\
    int x0 = *(stack0--);			\
    int y0 = *(stack0--);			\
    *(++stack0) = ((x0!=0) op (y0!=0));		\
    if(stack1){					\
      int x1 = *(stack1--);			\
      int y1 = *(stack1--);			\
      *(++stack1) = ((x1!=0) op (y1!=0));	\
    }						\
    break;}

#define ergm_LBINFUN(fun)			\
  {						\
    int x0 = *(stack0--);			\
    int y0 = *(stack0--);			\
    *(++stack0) = fun(x0!=0, y0!=0);		\
    if(stack1){					\
      int x1 = *(stack1--);			\
      int y1 = *(stack1--);			\
      *(++stack1) = fun(x1!=0, y1!=0);		\
    }						\
    break;}

static inline int ergm_LayerLogic2(Vertex ltail, Vertex lhead, // Dyad whose value/change to evaluate within the logical layer.
				   Vertex ttail, Vertex thead, // Dyad to toggle on LHS network.
				   StoreLayerLogic *ll, // Layer Logic
				   LayerLogicTask change
				  ){
  int *commands = ll->commands;
  // What gets looked up?
  Vertex lt = ltail, lh = lhead;
  // What gets toggled?
  Vertex tlt = ML_IO_TAIL(ll, ttail), tlh = ML_IO_HEAD(ll, thead), tl = ML_LID_TAIL(ll, ttail);
  // Is the dyad being toggled the same one as being looked up?
  unsigned int t_th = lt==tlt && lh==tlh, t_ht = ll->need_ht && lt==tlh && lh==tlt;

  int *stack0 = ll->stack0 - 1, // stack0 and stack1 always point to the top element (if any)
    *stack1 = change!=LL_ASIS && (t_th||t_ht)? ll->stack1 - 1 : NULL;  // Don't bother with stack1 if toggle can't affect focus dyad.

  int com;
  while((com = *(commands++)) != ML_STOP){
    switch(com){
    case 0:{
      int x0 = *(commands++);
      *(++stack0) = x0;
      if(stack1){
	*(++stack1) = x0;
      }
      break;}
    case -1:ergm_LUNOP(!)
    case -2:ergm_LBINOP(&&)
    case -3:ergm_LBINOP(||)
    case -4:ergm_LBINFUN(XOR)
    case -5:ergm_BINOP(==)
    case -6:ergm_BINOP(!=)
    case -7:ergm_BINOP(<)
    case -8:ergm_BINOP(>)
    case -9:ergm_BINOP(<=)
    case -10:ergm_BINOP(>=)
    case -11:ergm_BINOP(+)
    case -12:ergm_BINOP(-)
    case -13:ergm_BINOP(*)
    case -14:ergm_BINOP(/)
    case -15:ergm_BINOP(%)
    case -16:ergm_UNOP(-)
    case -17:ergm_UNFUN(abs)
    case -18:ergm_BINFUN(pow)
    case -19:ergm_BINFUN(fround)
    case -20:ergm_UNFUN(sign)
    case -21:{ // Reverse direction
      Vertex l = *(commands++);
      unsigned int x0;
      // If the physical layer is symmetrized, then only look at the upper triangle.
      if(ll->symm && ll->symm[l]!=0) x0 = ML_IGETWT(ll, l, MIN(lt,lh), MAX(lt,lh));
      else x0 = ML_IGETWT(ll, l, lh, lt);
      *(++stack0) = x0;
      if(stack1){
	unsigned int x1 = tl==l && (t_ht || (ll->symm && ll->symm[l]!=0 && t_th))? !(x0!=0) : x0;
	*(++stack1) = x1;
      }
      break;}
    default:{
      Vertex l = com; 
      unsigned int x0;
      // If the physical layer is symmetrized, then only look at the upper triangle.
      if(ll->symm && ll->symm[l]!=0) x0 = ML_IGETWT(ll, l, MIN(lt,lh), MAX(lt,lh));
      else x0 = ML_IGETWT(ll, l, lt, lh);
      *(++stack0) = x0;
      if(stack1){
	unsigned int x1 = tl==l && (t_th || (ll->symm && ll->symm[l]!=0 && t_ht))? !(x0!=0) : x0;
	*(++stack1) = x1;
      }
      break;}
    }
  }

  if(t_th||t_ht){
    switch(change){
    case LL_DIFF: return (int)(*stack1!=0) - (int)(*stack0!=0);
    case LL_ENCODE: return (*stack0!=0) | ((*stack1!=0)<<1);
    case LL_POST: return (*stack1!=0);
    default: return (*stack0!=0);
    }
  }else{
    switch(change){
    case LL_DIFF: return 0;
    case LL_ENCODE: return (*stack0!=0) | ((*stack0!=0)<<1);
    case LL_POST: return (*stack0!=0);
    default: return (*stack0!=0);
    }
  }
}

static inline int ergm_LayerLogic(Vertex tail, Vertex head, // Dyad to toggle and evaluate on LHS network.
				   StoreLayerLogic *ll, // Layer Logic
				   LayerLogicTask change
				  ){
  return ergm_LayerLogic2(ML_IO_TAIL(ll, tail), ML_IO_HEAD(ll, head), tail, head, ll, change);
}

// change = LL_ENCODE -> asis_th + toggled_th*2 + asis_ht*4 + toggled_ht*8
static inline unsigned int ergm_LayerLogic_affects(Vertex ttail, Vertex thead, // Dyad to toggle on LHS network.
						   StoreLayerLogic *ll, // Layer Logic
						   LayerLogicTask change,
						   Vertex *atails, Vertex *aheads){
  unsigned int nt = 0;
  Vertex lt = ML_IO_TAIL(ll, ttail), lh = ML_IO_HEAD(ll, thead);
  if(change==LL_ENCODE){
    return ergm_LayerLogic2(lt, lh, ttail, thead, ll, LL_ENCODE) | (ll->need_ht ? ergm_LayerLogic2(lt, lh, ttail, thead, ll, LL_ENCODE)<<2 : 0);
  }else{
    if(ergm_LayerLogic2(lt, lh, ttail, thead, ll, change)){
      if(atails) atails[nt] = lt;
      if(aheads) aheads[nt] = lh;
      nt++;
    }
    if(ll->need_ht && ergm_LayerLogic2(lh, lt, ttail, thead, ll, change)){
      if(atails) atails[nt] = lh;
      if(aheads) aheads[nt] = lt;
      nt++;
    }
  return nt;
  }
}
				  

#undef ergm_UNOP
#undef ergm_UNFUN
#undef ergm_BINOP
#undef ergm_BINFUN

#endif // _ERGM_MULTI_CHANGESTAT_MULTILAYER_H_
